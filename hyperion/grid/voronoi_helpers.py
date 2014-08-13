from astropy import log as logger


class voronoi_grid(object):
    '''
    A Voronoi grid.

    This class is used to compute properties of a 3D Voronoi grid useful
    in RT simulations, such as site neighbours lists, cell volumes, vertices
    bounding boxes, etc.

    The class needs at least two parameters for initialisation: the sites list
    and a domain representing a rectangular box containing the sites.

    Additional parameters include a flag to signal whether the vertices of the
    Voronoi cells should be returned or not, and two arguments defining
    additional (optional) walls within the domain.

    Internally, the Voronoi tessellation is computed via the voro++ library
    (http://math.lbl.gov/voro++/). The properties of the Voronoi grid
    can be retrieved via the ``neighbours_table`` class property.

    Parameters
    ----------
    sites : np.ndarray of floats
        The array of sites, must be a sequence of 3-dimensional vectors
        representing the position of each site.
    domain : np.ndarray of floats
        The domain of the grid, must be a 3x2 array corresponding to the
        min/max values in the three coordinates.
    with_vertices : boolean
        If ``True``, the vertices of the Voronoi cells will be computed.
        The vertices of the cells are not needed for RT simulations,
        and they add a considerable overhead in terms of memory requirements.
        Required for plotting and useful for debugging.
    wall : a string or ``None``
        If not ``None``, it must be one of
        ``['sphere','cylinder','cone','plane']``. It will define an additional
        wall of the desired shape within the domain.
    wall_args : a tuple of floats or ``None``
        Meaningful only if ``wall`` is not None, it represents the parameters
        to be used to construct the wall (e.g., in case of a spherical wall,
        the centre position and the radius of the sphere). Different construction
        arguments are required for different types of wall, see the voro++
        documentation for more information.
    '''

    def __init__(self, sites, domain, with_vertices=False, wall=None, wall_args=None, verbose=False):
        import numpy as np
        from ._voronoi_core import _voropp_wrapper
        from astropy.table import Table

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
        if sites.shape[1] != 3:
            raise ValueError(
                'the sites must be defined in 3-dimensional space')
        # Check that the domain has meaningful values.
        for limit in domain:
            if not limit[0] < limit[1]:
                raise ValueError('invalid domain')
        # Check that the input sites fall within the domain.
        for site in sites:
            for coord, limit in zip(site, domain):
                if coord <= limit[0] or coord >= limit[1]:
                    raise ValueError('a site is outside the domain')
        if not isinstance(with_vertices, bool):
            raise TypeError(
                'the \'with_vertices\' parameter must be a boolean')
        # Wall checks.
        allowed_walls = ['sphere','cylinder','cone','plane']
        if not wall is None and not wall in allowed_walls:
            raise ValueError('the \'wall\' parameter must be None or one of ' + str(allowed_walls))
        if not wall_args is None and (not isinstance(wall_args,tuple) or not all([isinstance(_,float) for _ in wall_args])):
            raise ValueError('the \'wall_args\' parameter must be None or a tuple of floats')

        # Redefine wall params in order to pass it to the C++ routine.
        wall = "" if wall is None else wall
        wall_args = () if wall_args is None else wall_args

        # Store the vertices flag.
        self._with_vertices = with_vertices

        logger.info("Computing the tessellation via voro++")
        tup = _voropp_wrapper(sites, domain, with_vertices, wall, wall_args, 1 if verbose else 0)
        if with_vertices:
            t = Table([sites, tup[0], tup[1], tup[2], tup[3], tup[4]],
                      names=('coordinates', 'neighbours', 'volume', 'bb_min', 'bb_max', 'vertices'))
        else:
            t = Table([sites, tup[0], tup[1], tup[2], tup[3]],
                      names=('coordinates', 'neighbours', 'volume', 'bb_min', 'bb_max'))
        self._neighbours_table = t

    # Getter for the neighbours/volume/bb table.
    @property
    def neighbours_table(self):
        from copy import deepcopy
        return deepcopy(self._neighbours_table)

    def plot(self):
        '''
        Plot a 3D visualisation of the Voronoi grid using mayavi.

        This method requires mayavi to be installed and also needs the vertices
        information to be available (see the class constructor).

        Note that in order for this method to work in an interactive IPython session,
        a series of environment variables and proper switches need to be used
        depending on your system configuration. For instance, on a Linux machine
        with PyQt4 and a recent IPython version, the following bash startup
        command for IPython can be used:
        ``ETS_TOOLKIT=qt4 QT_API=pyqt ipython --gui=qt4``
        This sets both the mayavi and the IPython GUI toolkit to qt4, and the ``QT_API``
        variable is used to specify that we want the ``pyqt`` API (as opposed to the
        ``pyside`` alternative API - PySide is an alternative implementation of PyQt).

        It should be possible to get this method working on different configurations,
        but the details will be highly system-dependent.
        '''

        if not self._with_vertices:
            raise ValueError(
                'the class must be constructed with \'with_vertices=True\' in order to support plotting')

        import numpy as np
        try:
            from tvtk.api import tvtk
            from mayavi.api import Engine
            from mayavi import mlab
            from mayavi.sources.vtk_data_source import VTKDataSource
            from mayavi.modules.surface import Surface
            from mayavi.modules.scalar_cut_plane import ScalarCutPlane
        except ImportError:
            raise ImportError(
                'the plot method requires Mayavi, please make sure it is correctly installed')

        # Shortcut.
        vertices = self._neighbours_table['vertices']

        # This is a list of all the vertices composing all voronoi cells.
        # points = [[x1,y1,z1],[x2,y2,z2],...]
        points = []
        # Array to describe each voronoi cell in terms of the points list above. E.g.,
        # cells = [4,0,1,2,3,5,4,5,6,7,8]
        # This describes two cells, the first with 4 vertices whose indices in the points array
        # are 0,1,2,3, the second with 5 vertices whose indices are 4,5,6,7,8.
        cells = []
        cur_cell_idx = 0
        # Indices in the cells array where each new cell starts. In the example above,
        # offset = [0,5]
        offset = []
        cur_offset = 0
        # Array of cell types. Cells will all be of the same type.
        cell_types = []

        # Build the above quantities.
        for v in vertices:
            # Drop the empty vertices coordinates, signalled by NaN.
            arr = v[~np.isnan(v)]
            assert(len(arr) % 3 == 0)
            tmp = np.split(arr, len(arr) / 3)
            # Append the vertices.
            points = points + tmp
            # Append the cell description.
            cells = cells + \
                [len(tmp)] + range(cur_cell_idx, cur_cell_idx + len(tmp))
            cur_cell_idx += len(tmp)
            # Append the offset info.
            offset.append(cur_offset)
            cur_offset += len(tmp) + 1
            # Append the cell type.
            cell_types.append(tvtk.ConvexPointSet().cell_type)

        # Cache the sites' positions.
        sites_arr = self._neighbours_table['coordinates']

        # Setup the Mayavi engine and figure.
        e = Engine()
        e.start()
        fig = mlab.figure(engine=e)

        # Plot the sites.
        mlab.points3d(
            sites_arr[:, 0], sites_arr[:, 1], sites_arr[:, 2], figure=fig)

        # Plot the cells with coloured surfaces.
        # This is just an array of scalars to assign a "temperature" to each cell vertex, which will be
        # used for coloring purposes.
        temperature = np.arange(0, len(points) * 10, 10, 'd')
        # Initialise the array of cells.
        cell_array = tvtk.CellArray()
        cell_array.set_cells(len(vertices), np.array(cells))
        # Initialise the unstructured grid object.
        ug = tvtk.UnstructuredGrid(points=np.array(points))
        ug.set_cells(np.array(cell_types), np.array(offset), cell_array)
        ug.point_data.scalars = temperature
        ug.point_data.scalars.name = 'temperature'
        # Create a data source from the unstructured grid object.
        src = VTKDataSource(data=ug)
        # Add the source to the engine.
        e.add_source(src)
        # Create a surface object with opacity 0.5
        surf = Surface()
        surf.actor.property.opacity = 0.5
        # Add the surface object to the engine.
        e.add_module(surf)
        # Add a cut plane as well.
        e.add_module(ScalarCutPlane())

        # Create another representation of the grid, this time using only white wireframe
        # to highlight to shape of the cells.
        # Rebuild the ug.
        ug = tvtk.UnstructuredGrid(points=np.array(points))
        ug.set_cells(np.array(cell_types), np.array(offset), cell_array)
        src = VTKDataSource(data=ug)
        e.add_source(src)
        surf = Surface()
        surf.actor.property.representation = 'wireframe'
        e.add_module(surf)
        cp = ScalarCutPlane()
        e.add_module(cp)
