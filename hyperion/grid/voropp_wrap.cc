#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "voro++/voro++.hh"

extern "C" const char * hyperion_voropp_wrap(int **sparse_neighbours, int **neigh_pos, int *nn, double **volumes, double **bb_min, double **bb_max, double **vertices,
                                             int *max_nv, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                                             double const *points, int npoints, int with_vertices, const char *wall_str, const double *wall_args_arr, int n_wall_args,
                                             int with_sampling, int n_samples, double **sample_points, int **sampling_idx, int *tot_samples, int min_cell_samples,
                                             int seed, int verbose);

using namespace voro;

// Number of average particles per block for good performance,
// determined experimentally.
static const double particle_block = 5.;

// Simple smart pointer that uses std::free to deallocate.
template <typename T>
class ptr_raii
{
    public:
        explicit ptr_raii(T *ptr):m_ptr(ptr) {}
        ~ptr_raii()
        {
            if (m_ptr) {
                std::free(m_ptr);
            }
        }
        T *release()
        {
            T *retval = m_ptr;
            m_ptr = 0;
            return retval;
        }
        T *get()
        {
            return m_ptr;
        }
    private:
        // Remove copy ctor and assignment operator.
        ptr_raii(const ptr_raii &);
        ptr_raii &operator=(const ptr_raii &);
    private:
        T *m_ptr;
};

// Functor to extract the max number of neighbours/vertices.
template <typename T>
static inline bool size_cmp(const std::vector<T> &a, const std::vector<T> &b)
{
    return a.size() < b.size();
}

// Global string used for reporting errors back to Python.
static std::string error_message;

// Compute the volume of a tetrahedron.
template <typename It>
static inline double tetra_volume(It v0, It v1, It v2, It v3)
{
    double mat[9];
    for (int i = 0; i < 3; ++i) {
        mat[i] = v1[i]-v0[i];
    }
    for (int i = 0; i < 3; ++i) {
        mat[3+i] = v2[i]-v0[i];
    }
    for (int i = 0; i < 3; ++i) {
        mat[6+i] = v3[i]-v0[i];
    }
    double a = mat[0], b = mat[1], c = mat[2];
    double d = mat[3], e = mat[4], f = mat[5];
    double g = mat[6], h = mat[7], i = mat[8];
    return std::abs(a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)) / 6.;
}

// http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html
// http://dl.acm.org/citation.cfm?id=378603
template <typename Ptr, typename It>
static inline void sample_point_in_tetra(Ptr res,It p0, It p1, It p2, It p3)
{
    // Three random numbers between zero and 1.
    double s = std::rand()/(RAND_MAX + 1.0);
    double t = std::rand()/(RAND_MAX + 1.0);
    double u = std::rand()/(RAND_MAX + 1.0);

    if (s + t > 1.0) {
        s = 1.0 - s;
        t = 1.0 - t;
    }

    if (t + u > 1.0) {
        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    } else if (s + t + u > 1.0) {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }

    const double a = 1.0 - s - t - u;

    for (int i = 0; i < 3; ++i) {
        res[i] = p0[i]*a+p1[i]*s+p2[i]*t+p3[i]*u;
    }
}

// Main wrapper called from cpython.
const char *hyperion_voropp_wrap(int **sparse_neighbours, int **neigh_pos, int *nn, double **volumes, double **bb_min, double **bb_max, double **vertices,
                                 int *max_nv, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double const *points,
                                 int nsites, int with_vertices, const char *wall_str, const double *wall_args_arr, int n_wall_args, int with_sampling, int n_samples,
                                 double **sample_points, int **sampling_idx, int *tot_samples, int min_cell_samples, int seed, int verbose)
{
    // We need to wrap everything in a try/catch block as exceptions cannot leak out to C.
    try {

    // Set the random seed.
    std::srand(static_cast<unsigned>(seed));

    // Volume of the domain.
    const double dom_vol = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);

    // Total number of blocks we want.
    const double nblocks = nsites / particle_block;

    // Average block edge.
    const double block_edge = cbrt(nblocks);

    // Average edge length of the domain.
    const double vol_edge = cbrt((xmax - xmin) * (ymax - ymin) * (zmax - zmin));

    // The number of grid blocks across each coordinate will be proportional
    // to the dimension of the domain in that coordinate. The +1 is to account for rounding
    // and to make sure that we always have at least 1 block.
    const int nx = (int)((xmax - xmin) / vol_edge * block_edge) + 1;
    const int ny = (int)((ymax - ymin) / vol_edge * block_edge) + 1;
    const int nz = (int)((zmax - zmin) / vol_edge * block_edge) + 1;

    if (verbose) {
        std::cout << "Number of sites: " << nsites << '\n';
        std::cout << "Domain: [" << xmin << ',' << xmax << "] [" << ymin << ',' << ymax << "] [" << zmin << ',' << zmax << "]\n";
        std::cout << "Initialising with the following block grid: " << nx << ',' << ny << ',' << nz << '\n';
        std::cout << std::boolalpha;
        std::cout << "Vertices: " << bool(with_vertices) << '\n';
        std::cout << "With sampling: " << bool(with_sampling) << '\n';
    }

    // Prepare the output quantities.
    // Neighbour list.
    std::vector<std::vector<int> > n_list(nsites);
    // List of vertices.
    std::vector<std::vector<double> > vertices_list;
    if (with_vertices) {
        vertices_list.resize(nsites);
    }
    // Volumes.
    ptr_raii<double> vols(static_cast<double *>(std::malloc(sizeof(double) * nsites)));
    // Bounding boxes.
    ptr_raii<double> bb_m(static_cast<double *>(std::malloc(sizeof(double) * nsites * 3)));
    ptr_raii<double> bb_M(static_cast<double *>(std::malloc(sizeof(double) * nsites * 3)));

    // Initialise the voro++ container.
    container con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
                  false,false,false,8);
    for(int i = 0; i < nsites; ++i) {
            con.put(i,points[i*3],points[i*3 + 1],points[i*3 + 2]);
    }

    // Initialise the looping variables and the temporary cell object used for computation.
    voronoicell_neighbor c;
    c_loop_all vl(con);
    int idx;
    double tmp_min[3],tmp_max[3];
    // Vector to store temporarily the list of vertices coordinates triplets. Used only if vertices
    // are not requested (otherwise, the vertices are stored directly in the output array).
    std::vector<double> tmp_v;
    // Site position and radius (r is unused).
    double x,y,z,r;
    // These quantities are sampling-related.
    // List of faces for each cell. Format described here:
    // http://math.lbl.gov/voro++/examples/polygons/
    std::vector<int> f_vert;
    // List of tetrahedra vertices indices: 4 elements per tet.
    std::vector<int> t_vert;
    // Cumulative volumes of the tetrahedra.
    std::vector<double> c_vol;
    // Vector of vectors of sampling points for each cell.
    std::vector<std::vector<double> > vs_points;
    vs_points.resize(nsites);
    // Vectors of the indices of the sampling points.
    std::vector<int> vs_points_idx;
    vs_points_idx.push_back(0);
    // Init the total number of samples.
    *tot_samples = 0;
    // A vector of random values in the [0,1[ range used for the part of the sampling proportional
    // to the cell volumes.
    std::vector<double> random_values;
    for (int i = 0; i < n_samples; ++i) {
        random_values.push_back(std::rand()/(RAND_MAX + 1.0));
    }
    // Sort it in ascending order.
    std::sort(random_values.begin(),random_values.end());
    // Cumulative volume vector for the cells.
    std::vector<double> cc_vol;
    cc_vol.push_back(0);

    // Init the total number of neighbours.
    *nn = 0;

    // Loop over all particles and compute the desired quantities.
    if(vl.start()) {
        do {
            // Get the id and position of the site being considered.
            vl.pos(idx,x,y,z,r);
            std::vector<double> *tmp_vertices = with_vertices ? &(vertices_list[idx]) : &tmp_v;
            // Compute the voronoi cell.
            con.compute_cell(c,vl);
            // Compute the neighbours.
            c.neighbors(n_list[idx]);
            // Update the number of neighbours.
            *nn += n_list[idx].size();
            // Volume.
            vols.get()[idx] = c.volume();
            // Compute bounding box. Start by asking for the vertices.
            c.vertices(x,y,z,*tmp_vertices);
            // Init min/max bb.
            std::copy(tmp_vertices->begin(),tmp_vertices->begin() + 3,tmp_min);
            std::copy(tmp_vertices->begin(),tmp_vertices->begin() + 3,tmp_max);
            for (unsigned long i = 1u; i < tmp_vertices->size() / 3u; ++i) {
                for (unsigned j = 0; j < 3; ++j) {
                    if ((*tmp_vertices)[i * 3 + j] < tmp_min[j]) {
                        tmp_min[j] = (*tmp_vertices)[i * 3 + j];
                    }
                    if ((*tmp_vertices)[i * 3 + j] > tmp_max[j]) {
                        tmp_max[j] = (*tmp_vertices)[i * 3 + j];
                    }
                }
            }
            // Copy the bounding box into the output array.
            std::copy(tmp_min,tmp_min + 3,bb_m.get() + idx * 3);
            std::copy(tmp_max,tmp_max + 3,bb_M.get() + idx * 3);
            // Computation of the sampling array, only if requested.
            if (!with_sampling) {
                continue;
            }
            // Clear tmp variables.
            t_vert.clear();
            c_vol.clear();
            double vol = 0;
            // Update the cumulative volume vector.
            cc_vol.push_back(cc_vol.back() + vols.get()[idx] / dom_vol);
            // Compute the faces of the cell.
            c.face_vertices(f_vert);
            int j = 0;
            while (j < f_vert.size()) {
                // Number of vertices for each face.
                const int nfv = f_vert[j];
                // We need to establish if the vertex with index 0 in the cell is
                // part of the current face. If that is the case, we skip
                // the face.
                if (std::find(f_vert.begin() + j + 1,f_vert.begin() + j + 1 + nfv,0) != f_vert.begin() + j + 1 + nfv) {
                    // Don't forget to update the counter...
                    j += nfv + 1;
                    continue;
                }
                // Now we need to build the list of tetrahedra vertices. The procedure:
                // - the first vertex is always the 0 vertex of the cell,
                // - the second vertex is always the first vertex of the current face.
                // - the other two vertices are taken as the k-th and k+1-th vertices of
                //   the face.
                for (int k = j + 2; k < j + nfv; ++k) {
                    t_vert.push_back(0);
                    t_vert.push_back(f_vert[j+1]);
                    t_vert.push_back(f_vert[k]);
                    t_vert.push_back(f_vert[k+1]);
                    // Volume of the tetrahedron.
                    double t_vol = tetra_volume(tmp_vertices->begin(),tmp_vertices->begin() + f_vert[j+1]*3,tmp_vertices->begin() + f_vert[k]*3,
                                                tmp_vertices->begin() + f_vert[k+1]*3);
                    // Update the cumulative volume and add it.
                    vol += t_vol;
                    c_vol.push_back(vol);
                }
                // Update the counter.
                j += nfv + 1;
            }
            // Now we need to establish how many points we need to randomly pick in this cell. The value must be
            // at least min_cell_samples and proportional to the cell volume. In order to do this, we need
            // to compute how many random samples in random_values fall within the current and previous cumulative
            // cell volume.
            typedef std::vector<double>::iterator it_type;
            // First element in random_values that is greater than or equal to the lower bound of the volume range.
            it_type lower = std::lower_bound(random_values.begin(),random_values.end(),*(cc_vol.end() - 2));
            // First element in random_values that is greater than the upper bound of the volume range.
            it_type upper = std::upper_bound(random_values.begin(),random_values.end(),*(cc_vol.end() - 1));
            // The number of elements to sample according to the cell volume will be the distance between upper and lower.
            int nc_samples = std::distance(lower,upper);
            // If less than min samples, overwrite it.
            if (nc_samples < min_cell_samples) {
                nc_samples = min_cell_samples;
            }
            // Now we need to select randomly tetras (with a probability proportional to their volume) and sample
            // uniformly inside them.
            const double c_factor = c_vol.back()/(RAND_MAX + 1.0);
            // Go with the sampling.
            for (int i = 0; i < nc_samples;) {
                const double r_vol = std::rand()*c_factor;
                std::vector<double>::iterator it = std::upper_bound(c_vol.begin(),c_vol.end(),r_vol);
                // It might be possible due to floating point madness that this actually goes past the end,
                // in such a case we just repeat the sampling.
                if (it == c_vol.end()) {
                    continue;
                }
                // Make space for the new point.
                vs_points[idx].push_back(0);
                vs_points[idx].push_back(0);
                vs_points[idx].push_back(0);
                // Get the index of the tetra that was selected.
                int t_idx = std::distance(c_vol.begin(),it);
                // Now we can go sample inside the selected tetrahedron.
                std::vector<int>::iterator t_it = t_vert.begin() + t_idx*4;
                sample_point_in_tetra(&vs_points[idx].back() - 2,tmp_vertices->begin()+(*t_it)*3,tmp_vertices->begin()+(*(t_it + 1))*3,
                                      tmp_vertices->begin()+(*(t_it + 2))*3,tmp_vertices->begin()+(*(t_it + 3))*3);
                ++i;
                // Update the total number of samples.
                ++(*tot_samples);
            }
        } while(vl.inc());
    }

    // Copy the sampling points in the final array, if requested.
    ptr_raii<double> spoints(with_sampling ? static_cast<double *>(std::malloc(sizeof(double) * (*tot_samples) * 3)) : NULL);
    ptr_raii<int> spoints_idx(with_sampling ? static_cast<int *>(std::malloc(sizeof(int) * (nsites + 1))) : NULL);
    if (with_sampling) {
        double *cur_sample = spoints.get();
        *(spoints_idx.get()) = 0;
        int acc = 0;
        for (int i = 0; i < nsites; ++i) {
            std::copy(vs_points[i].begin(),vs_points[i].end(),cur_sample);
            cur_sample += vs_points[i].size();
            acc += vs_points[i].size() / 3;
            *(spoints_idx.get() + i + 1) = acc;
        }
        // Free the vector.
        vs_points = std::vector<std::vector<double> >();
    }

    // The voro++ doc say that in case of numerical errors the neighbours list might not be symmetric,
    // that is, if 'a' is a neighbour of 'b' then 'b' might not be a neighbour of 'a'. We check and fix this
    // in the loop below.
    for (idx = 0; idx < nsites; ++idx) {
        for (unsigned j = 0u; j < n_list[idx].size(); ++j) {
            // Check only non-wall neighbours.
            if (n_list[idx][j] >= 0 && std::find(n_list[n_list[idx][j]].begin(),n_list[n_list[idx][j]].end(),idx) == n_list[n_list[idx][j]].end()) {
                n_list[n_list[idx][j]].push_back(idx);
            }
        }
    }

    // Allocate space for the sparse version of the neighbour indices.
    ptr_raii<int> sparse_neighs(static_cast<int *>(std::malloc(sizeof(int) * (*nn))));
    ptr_raii<int> sparse_neighs_indices(static_cast<int *>(std::malloc(sizeof(int) * (nsites + 1))));
    // Build the sparse version of the neighbours list.
    int cur_idx = 0;
    for (std::size_t i = 0u; i < n_list.size(); ++i) {
        const std::vector<int> &v = n_list[i];
        std::copy(v.begin(),v.end(),sparse_neighs.get() + cur_idx);
        *(sparse_neighs_indices.get() + i) = cur_idx;
        cur_idx += v.size();
    }
    // Last index.
    *(sparse_neighs_indices.get() + nsites) = cur_idx;

    if (with_vertices) {
        // Compute the max number of vertices coordinates.
        *max_nv = std::max_element(vertices_list.begin(),vertices_list.end(),size_cmp<double>)->size();
        if (verbose) std::cout << "Max number of vertices coordinates is: " << *max_nv << '\n';

        // Allocate space for flat array of vertices.
        ptr_raii<double> verts(static_cast<double *>(std::malloc(sizeof(double) * nsites * (*max_nv))));
        // Fill it in.
        for (idx = 0; idx < nsites; ++idx) {
            double *ptr = verts.get() + (*max_nv) * idx;
            std::copy(vertices_list[idx].begin(),vertices_list[idx].end(),ptr);
            // Fill empty elements with nan.
            std::fill(ptr + vertices_list[idx].size(),ptr + (*max_nv),std::numeric_limits<double>::quiet_NaN());
        }

        // Assign the output quantity.
        *vertices = verts.release();
    } else {
        *max_nv = 0;
    }

    // Assign the output quantities.
    *volumes = vols.release();
    *bb_min = bb_m.release();
    *bb_max = bb_M.release();
    *sample_points = spoints.release();
    *sampling_idx = spoints_idx.release();
    *sparse_neighbours = sparse_neighs.release();
    *neigh_pos = sparse_neighs_indices.release();

    return NULL;

    } catch (const std::exception &e) {
        error_message = std::string("A C++ exception was raised while calling the voro++ wrapper. The full error message is: \"") + e.what() + "\".";
        return error_message.c_str();
    }
}
