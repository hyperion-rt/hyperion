import multiprocessing as _mp

# Function used to disable ctrl+c in sub processes.


def _do_nothing(signum, frame):
    pass

# Parallel functor for the computation of the initial
# tessellations.


def _par_tess(t):
    import signal

    signal.signal(signal.SIGINT, _do_nothing)

    retval = t[0](t[1])

    return retval

# Parallel functor for the computation of the volumes
# and the bounding boxes.


def _par_vol_bb(t):
    import signal
    import numpy as np
    from scipy.spatial import Delaunay

    signal.signal(signal.SIGINT, _do_nothing)

    # Unpack the arguments from the input tuple.
    chunk, points, point_region, protruding_cells, regions, v_vertices, domain, idx = t
    start, end = chunk

    try:
        # Try to pin the process to a processor using the (optional) psutil module.
        # If anything fails, do nothing.
        import psutil
        import os

        p = psutil.Process(os.getpid())
        p.set_cpu_affinity([idx])
    except:
        pass

    ndim = len(domain)

    # Function to convert from site indices to region indices.
    sidx_to_ridx = lambda sidx: point_region[sidx]

    # Init the array of bounding boxes with zeroes.
    bb_arr = np.zeros([end - start, 2 * ndim])

    # Init the volumes array with -1s.
    vol_arr = np.empty([end - start], dtype=np.dtype(float))
    vol_arr.fill(-1.)

    # Main loop for the computation of the volumes and bounding boxes.
    for i in range(start, end):
        ridx = sidx_to_ridx(i)

        # Protruding cells will have a volume of -1 and bounding box
        # of [[0 ... 0],[0 ... 0]].
        if not ridx in protruding_cells:
            cell = regions[ridx]
            vertices = v_vertices[cell]

            # Decompose the Voronoi region into simplices.
            dtess = Delaunay(vertices)

            # Accumulate the volumes of the simplices.
            volume = 0.
            for simplex in dtess.simplices:
                volume += _simplex_volume(dtess.points[simplex])
            vol_arr[i - start] = volume

            # Bounding box.
            # Initialise min/max with the coordinates of the first vertex.
            bb_arr[i - start][0:ndim] = vertices[0]
            bb_arr[i - start][ndim:] = vertices[0]

            # Iterate on the remaining vertices and update the minmax values as
            # needed.
            for vertex in vertices[1:]:
                for j in range(0, len(vertex)):
                    if vertex[j] < bb_arr[i - start][j]:
                        bb_arr[i - start][j] = vertex[j]
                    if vertex[j] > bb_arr[i - start][ndim + j]:
                        bb_arr[i - start][ndim + j] = vertex[j]

    return vol_arr, bb_arr

# Compute the volume of a simplex.
# http://en.wikipedia.org/wiki/Simplex#Volume


def _simplex_volume(simplex):
    n = len(simplex) - 1

    # Special casing for 3d, use the C version.
    if n == 3:
        from ._voronoi_core import _simplex3d_volume
        return _simplex3d_volume(simplex)

    from math import gamma
    import numpy as np

    matrix = np.empty((n, n))
    for i in range(0, n):
        matrix[i] = simplex[i + 1] - simplex[0]

    det = np.linalg.det(matrix)
    return abs(det / gamma(n + 1))


class voronoi_grid(object):

    def __init__(self, sites, domain, ncpus=_mp.cpu_count()):
        from scipy.spatial import Delaunay, Voronoi
        import numpy as np
        from copy import deepcopy

        # Validate input.
        if not isinstance(sites, np.ndarray) or sites.dtype.kind != 'f':
            raise TypeError('input sites list must be a NumPy array of floats')
        if not isinstance(domain, np.ndarray) or domain.dtype.kind != 'f':
            raise TypeError('input domain must be a NumPy array of floats')
        if len(sites.shape) != 2:
            raise ValueError('the input sites list must be a 2d array')
        if len(domain.shape) != 2:
            raise ValueError('the input domain must be a 2d array')
        if sites.shape[1] != domain.shape[0]:
            raise ValueError(
                'the size of the domain must match the dimensionality of the sites')
        if domain.shape[1] != 2:
            raise ValueError(
                'the domain must be an array of [min,max] elements')
        if sites.shape[1] < 2:
            raise ValueError(
                'the sites must be defined in a space at least 2-dimensional')
        # Check that the domain has meaningful values.
        for limit in domain:
            if not limit[0] < limit[1]:
                raise ValueError('invalid domain')
        # Check that the input sites fall within the domain.
        for site in sites:
            for coord, limit in zip(site, domain):
                if coord < limit[0] or coord > limit[1]:
                    raise ValueError('a site is outside the domain')

        # Setup the process pool.
        self._pool = _mp.Pool(ncpus)
        self._ncpus = ncpus

        # Store the domain.
        self._domain = deepcopy(domain)

        # Build the two tessellations.
        # NOTE: the two tessellations use different types of indexing:
        # - the Delaunay uses the indices of the input sites,
        # - the Voronoi refers to the indices of the regions created from the sites,
        #   which will be in general unrelated to the indices of the input sites.
        # We provide below a couple of utility functions to convert between the two indexing schemes.
        # Note that the user of the class cares about the indices of the input sites, so everything
        # user-facing should be presented in that convention.
        self._del_tess, self._vor_tess = self._pool.map(
            _par_tess, [(Delaunay, sites), (Voronoi, sites)])

        # Close the pool for now, free resources.
        self._pool.close()
        self._pool.join()

        # Build the dict for the mapping from region idx to site idx.
        self._r_to_s_dict = dict(
            [(t[1], t[0]) for t in enumerate(self._vor_tess.point_region)])

        # Compute the initial list of neighbours. This is stored in sites
        # indices.
        self._nl = self._compute_neighbours()

        # Identify the cells that are protruding from the box (including infinite cells).
        # Here we are using the Voronoi regions, so everything is indexed
        # according to the regions.
        self._protruding_cells = self._compute_protruding_cells()

        # Compute and store the table.
        self._neighbours_table = self._compute_neighbours_table()

        # NOTE: this is used for the debug plotting methods, disabled for now.
        # The wall cells will be all the neighbours of the protruding ones, excluding the
        # protruding ones themselves.
        #wall_cells = set()
        # for cell in protruding_cells:
        # for neighbour in self._nl[self._ridx_to_sidx(cell)]:
        # if neighbour != -1 and not self._sidx_to_ridx(neighbour) in protruding_cells:
        # wall_cells.add(self._sidx_to_ridx(neighbour))
        #self._new_regions = regions_copy
        #self._wall_cells = wall_cells

    # Convert site idx (i.e., indexing with respect to the input sites, used in the Delaunay tessellation
    # and in the neighbours list) to the indexing of Voronoi regions as
    # resulting from the Voronoi tessellation.
    def _sidx_to_ridx(self, sidx):
        return self._vor_tess.point_region[sidx]

    # Opposite of the above.
    def _ridx_to_sidx(self, ridx):
        return self._r_to_s_dict[ridx]

    # Computation of protruding cells.
    def _compute_protruding_cells(self):
        #from copy import deepcopy
        from ._voronoi_core import _region_in_domain
        protruding_cells = set()
        #regions_copy = deepcopy(self._vor_tess.regions)
        for i in range(len(self._vor_tess.regions)):
            region = self._vor_tess.regions[i]
            # if -1 in region or not self._region_in_domain(region):
            if -1 in region or not _region_in_domain(region, self._vor_tess.vertices, self._domain):
                #regions_copy[i] = []
                protruding_cells.add(i)
        return protruding_cells

    # Compute neighbours list.
    def _compute_neighbours(self):
        import numpy as np
        from ._voronoi_core import _neighbours_list_loop
        neigh_list = [set() for i in range(0, len(self._del_tess.points))]

        # Run the main loop from C.
        _neighbours_list_loop(self._del_tess.simplices, neigh_list)

        # Transform the list into a NumPy array.
        n, m = len(neigh_list), len(max(neigh_list, key=lambda l: len(l)))
        # Fill with -1 values.
        retval = np.empty([n, m], dtype=np.int32)
        retval.fill(-1)
        for i in range(n):
            retval[i][:len(neigh_list[i])] = list(neigh_list[i])

        return retval

    # A couple of plotting methods used for debug.
    def _plot_original(self):
        from scipy.spatial import voronoi_plot_2d
        voronoi_plot_2d(self._vor_tess)

    def _plot_new(self):
        from matplotlib.pylab import plot, figure
        import numpy as np

        if len(self._domain) != 2:
            raise ValueError('only 2-dimensional plotting is implemented')

        #fig = figure()
        plotted_lines = set()

        # Plot a single region.
        def plot_region(idx):
            region = self._new_regions[idx]
            # Do nothing for empty regions.
            if region == []:
                return

            assert(len(region) >= 3)

            cur_vertex = region[0]
            for vertex_idx in region[1:] + [cur_vertex]:
                line = tuple(sorted([cur_vertex, vertex_idx]))
                if not line in plotted_lines or idx in self._wall_cells:
                    plotted_lines.add(line)
                    arr = np.array(
                        [self._vor_tess.vertices[line[0]], self._vor_tess.vertices[line[1]]]).transpose()
                    if idx in self._wall_cells:
                        plot(arr[0], arr[1], 'r')
                    else:
                        plot(arr[0], arr[1], 'b')
                cur_vertex = vertex_idx

        # Plot all the new regions.
        for i in range(len(self._new_regions)):
            plot_region(i)

    # Compute and return the neighbours/volume/bb table.
    def _compute_neighbours_table(self):
        from astropy.table import Table
        import numpy as np
        from scipy.spatial import Delaunay

        # Revive the pool.
        self._pool = _mp.Pool(self._ncpus)

        ndim = len(self._domain)

        ntot = len(self._del_tess.points)
        chunk_size = ntot // self._ncpus
        chunks = []
        for i in range(0, self._ncpus - 1):
            chunks.append((i * chunk_size, (i + 1) * chunk_size))
        # Last chunk.
        chunks.append(((self._ncpus - 1) * chunk_size, ntot))

        vblists = self._pool.map(_par_vol_bb, [(chunk, self._vor_tess.points, self._vor_tess.point_region,
                                                self._protruding_cells, self._vor_tess.regions, self._vor_tess.vertices, self._domain, idx) for chunk,idx in zip(chunks,range(self._ncpus))])

        self._pool.close()
        self._pool.join()

        vlist = [t[0] for t in vblists]
        blist = [t[1] for t in vblists]

        vol_arr = np.hstack(tuple(vlist))

        # Bounding boxes.
        bb_arr = np.vstack(tuple(blist))

        # Build the table, including the sites coordinates.
        t = Table([self._vor_tess.points, self._nl, vol_arr, bb_arr[:, 0:ndim], bb_arr[
                  :, ndim:]], names=('coordinates', 'neighbours', 'volume', 'bb_min', 'bb_max'))

        self._pool.close()
        self._pool.join()

        return t

    # Getter for the neighbours/volume/bb table.
    @property
    def neighbours_table(self):
        from copy import deepcopy
        return deepcopy(self._neighbours_table)
