class voronoi_grid(object):
    def __init__(self,sites,domain):
        from scipy.spatial import Delaunay, Voronoi
        import numpy as np
        from copy import deepcopy
        # Validate input.
        if not isinstance(sites,np.ndarray) or sites.dtype.kind != 'f':
            raise TypeError('input sites list must be a NumPy array of floats')
        if not isinstance(domain,np.ndarray) or domain.dtype.kind != 'f':
            raise TypeError('input domain must be a NumPy array of floats')
        if len(sites.shape) != 2:
            raise ValueError('the input sites list must be a 2d array')
        if len(domain.shape) != 2:
            raise ValueError('the input domain must be a 2d array')
        if sites.shape[1] != domain.shape[0]:
            raise ValueError('the size of the domain must match the dimensionality of the sites')
        if domain.shape[1] != 2:
            raise ValueError('the domain must be an array of [min,max] elements')
        if sites.shape[1] < 2:
            raise ValueError('the sites must be defined in a space at least 2-dimensional')
        # Check that the domain has meaningful values.
        for limit in domain:
            if not limit[0] < limit[1]:
                raise ValueError('invalid domain')
        # Check that the input sites fall within the domain.
        for site in sites:
            for coord,limit in zip(site,domain):
                if coord < limit[0] or coord > limit[1]:
                    raise ValueError('a site is outside the domain')
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
        self._del_tess = Delaunay(sites)
        self._vor_tess = Voronoi(sites)
        # Build the dict for the mapping from region idx to site idx.
        self._r_to_s_dict = dict([(t[1],t[0]) for t in enumerate(self._vor_tess.point_region)])
        # Compute the initial list of neighbours. This is stored in sites indices.
        self._nl = self._compute_neighbours()
        # Identify the cells that are protruding from the box (including infinite cells).
        # Here we are using the Voronoi regions, so everything is indexed according to the regions.
        protruding_cells = set()
        regions_copy = deepcopy(self._vor_tess.regions)
        for i in range(len(self._vor_tess.regions)):
            region = self._vor_tess.regions[i]
            if -1 in region or not self._region_in_domain(region):
                regions_copy[i] = []
                protruding_cells.add(i)
        self._protruding_cells = protruding_cells
        # The wall cells will be all the neighbours of the protruding ones, excluding the
        # protruding ones themselves.
        wall_cells = set()
        for cell in protruding_cells:
            for neighbour in self._nl[self._ridx_to_sidx(cell)]:
                if not self._sidx_to_ridx(neighbour) in protruding_cells:
                    wall_cells.add(self._sidx_to_ridx(neighbour))
        self._new_regions = regions_copy
        self._wall_cells = wall_cells
    # Convert site idx (i.e., indexing with respect to the input sites, used in the Delaunay tessellation
    # and in the neighbours list) to the indexing of Voronoi regions as resulting from the Voronoi tessellation.
    def _sidx_to_ridx(self,sidx):
        return self._vor_tess.point_region[sidx]
    # Opposite of the above.
    def _ridx_to_sidx(self,ridx):
        return self._r_to_s_dict[ridx]
    def _compute_neighbours(self):
        neigh_list = [set() for i in range(0,len(self._del_tess.points))]
        for simplex in self._del_tess.simplices:
            for site_idx in simplex:
                for neigh_idx in simplex:
                    if neigh_idx != site_idx:
                        neigh_list[site_idx].add(neigh_idx)
        return neigh_list
    # Test if a Voronoi region is entirely within the domain.
    def _region_in_domain(self,region):
        for vertex_idx in region:
            vertex = self._vor_tess.vertices[vertex_idx]
            for coord,limit in zip(vertex,self._domain):
                if coord < limit[0] or coord > limit[1]:
                    return False
        return True
    def _plot_original(self):
        from scipy.spatial import voronoi_plot_2d
        voronoi_plot_2d(self._vor_tess)
    def _plot_new(self):
        if len(self._domain) != 2:
            raise ValueError('only 2-dimensional plotting is implemented')
        from matplotlib.pylab import plot, figure
        import numpy as np
        #fig = figure()
        plotted_lines = set()
        def plot_region(idx):
            region = self._new_regions[idx]
            # Do nothing for empty regions.
            if region == []:
                return
            assert(len(region) >= 3)
            cur_vertex = region[0]
            for vertex_idx in region[1:] + [cur_vertex]:
                line = tuple(sorted([cur_vertex,vertex_idx]))
                if not line in plotted_lines or idx in self._wall_cells:
                    plotted_lines.add(line)
                    arr = np.array([self._vor_tess.vertices[line[0]],self._vor_tess.vertices[line[1]]]).transpose()
                    if idx in self._wall_cells:
                        plot(arr[0],arr[1],'r')
                    else:
                        plot(arr[0],arr[1],'b')
                cur_vertex = vertex_idx
        for i in range(len(self._new_regions)):
            plot_region(i)
    # Compute the volume of a simplex.
    def _simplex_volume(self,simplex):
        from math import gamma
        import numpy as np
        n = len(simplex) - 1
        matrix = np.zeros((n,n))
        for i in range(0,n):
            matrix[i] = simplex[i + 1] - simplex[0]
        det = np.linalg.det(matrix)
        return abs(det / gamma(n + 1))
    def neighbours_table(self):
        from astropy.table import Table
        import numpy as np
        from scipy.spatial import Delaunay
        # Establish the maximum number of neighbours.
        max_nn = len(max(self._nl,key = lambda l: len(l)))
        # Create the array of neighbours. Padding will be indicated by the value "-1".
        n_array = np.array([[-1] * max_nn] * len(self._vor_tess.points))
        # Fill in the array of neighbours.
        for i in range(len(self._nl)):
            tmp_list = list(self._nl[i])
            n_array[i][0:len(tmp_list)] = tmp_list
        # Now onto the volumes. Protruding cells will have a volume of -1.
        vol_list = []
        for i in range(len(self._vor_tess.points)):
            if self._sidx_to_ridx(i) in self._protruding_cells:
                vol_list.append(-1.)
            else:
                cell = self._vor_tess.regions[self._sidx_to_ridx(i)]
                vertices = self._vor_tess.vertices[cell]
                dtess = Delaunay(vertices)
                volume = 0.
                for simplex in dtess.simplices:
                    volume += self._simplex_volume(dtess.points[simplex])
                vol_list.append(volume)
        # Bounding boxes.
        bb_list = []
        for i in range(len(self._vor_tess.points)):
            ridx = self._sidx_to_ridx(i)
            if ridx in self._protruding_cells:
                bb_list.append(np.zeros((3,2)))
            else:
                cell = self._vor_tess.regions[ridx]
                vertices = self._vor_tess.vertices[cell]
                bb = [np.array([coord,coord]) for coord in vertices[0]]
                for vertex in vertices[1:]:
                    for j in range(0,len(vertex)):
                        if vertex[j] < bb[j][0]:
                            bb[j][0] = vertex[j]
                        if vertex[j] > bb[j][1]:
                            bb[j][1] = vertex[j]
                bb_list.append(np.array(bb))
        # Sites coordinates.
        t = Table([self._vor_tess.points,n_array,np.array(vol_list),np.array(bb_list)],names=('coordinates','neighbours','volume','bounding_box'))
        return t
