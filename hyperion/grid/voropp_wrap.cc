#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "voro++/voro++.hh"

extern "C" const char * hyperion_voropp_wrap(int **neighbours, int *max_nn, double **volumes, double **bb_min, double **bb_max,
                  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                  double const *points, int npoints);

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

// Functor to extract the max number of neighbours.
static inline bool nlist_cmp(const std::vector<int> &a, const std::vector<int> &b)
{
    return a.size() < b.size();
}

// Global string used for reporting errors back to Python.
static std::string error_message;

const char *hyperion_voropp_wrap(int **neighbours, int *max_nn, double **volumes, double **bb_min, double **bb_max,
       double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double const *points,
       int nsites)
{
    // We need to wrap everything in a try/catch block as exceptions cannot leak out to C.
    try {

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

    std::cout << "Number of sites: " << nsites << '\n';
    std::cout << "Domain: [" << xmin << ',' << xmax << "] [" << ymin << ',' << ymax << "] [" << zmin << ',' << zmax << "]\n";
    std::cout << "Initialising with the following block grid: " << nx << ',' << ny << ',' << nz << '\n';

    // Prepare the output quantities.
    // Neighbour list.
    std::vector<std::vector<int> > n_list(nsites);
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
    int idx = 0;
    // Temporary vector to store vertices information.
    std::vector<double> tmp_vertices;
    double tmp_min[3],tmp_max[3];

    // Loop over all particles and compute the desired quantities.
    if(vl.start()) {
        do {
            // Compute the voronoi cell.
            con.compute_cell(c,vl);
            // Compute the neighbours.
            c.neighbors(n_list[idx]);
            // Volume.
            vols.get()[idx] = c.volume();
            // Compute bounding box. Start by asking for the vertices.
            c.vertices(vl.x(),vl.y(),vl.z(),tmp_vertices);
            // Init min/max bb.
            std::copy(tmp_vertices.begin(),tmp_vertices.begin() + 3,tmp_min);
            std::copy(tmp_vertices.begin(),tmp_vertices.begin() + 3,tmp_max);
            for (unsigned long i = 0; i < tmp_vertices.size() / 3u; ++i) {
                for (unsigned j = 0; j < 3; ++j) {
                    if (tmp_vertices[i * 3 + j] < tmp_min[j]) {
                        tmp_min[j] = tmp_vertices[i * 3 + j];
                    }
                    if (tmp_vertices[i * 3 + j] > tmp_max[j]) {
                        tmp_max[j] = tmp_vertices[i * 3 + j];
                    }
                }
            }
            // Copy the bounding box into the output array.
            std::copy(tmp_min,tmp_min + 3,bb_m.get() + idx * 3);
            std::copy(tmp_max,tmp_max + 3,bb_M.get() + idx * 3);
            ++idx;
        } while(vl.inc());   
    }

    // Compute the max number of neighbours.
    *max_nn = std::max_element(n_list.begin(),n_list.end(),nlist_cmp)->size();
    std::cout << "Max number of neighbours is: " << *max_nn << '\n';

    // Allocate space for flat array of neighbours.
    ptr_raii<int> neighs(static_cast<int *>(std::malloc(sizeof(int) * nsites * (*max_nn))));
    // Fill it in.
    for (idx = 0; idx < nsites; ++idx) {
        int *ptr = neighs.get() + (*max_nn) * idx;
        std::copy(n_list[idx].begin(),n_list[idx].end(),ptr);
        // Fill empty elements with INT_MAX.
        std::fill(ptr + n_list[idx].size(),ptr + (*max_nn),INT_MAX);
    }

    // Assign the output quantities.
    *volumes = vols.release();
    *bb_min = bb_m.release();
    *bb_max = bb_M.release();
    *neighbours = neighs.release();

    return NULL;

    } catch (const std::exception &e) {
        error_message = std::string("A C++ exception was raised while calling the voro++ wrapper. The full error message is: \"") + e.what() + "\".";
        return error_message.c_str();
    }
}
