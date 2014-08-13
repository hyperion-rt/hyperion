from __future__ import print_function, division

import struct
import hashlib
from copy import deepcopy

import h5py
import numpy as np

from astropy import log as logger
from astropy.table import Table

from ..util.meshgrid import meshgrid_nd
from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing, link_or_copy
from .grid_helpers import single_grid_dims


class VoronoiGrid(FreezableClass):
    '''
    A voronoi mesh.

    The mesh can be initialized by passing the x, y, and z coordinates of the
    points used to compute the mesh::

        >>> grid = Voronoi(x, y, z)

    where ``x``, ``y``, and ``z`` are 1-d sequences of point positions.

    :class:`~hyperion.grid.VoronoiGrid` objects may contain multiple
    quantities (e.g. density, specific energy). To access these, you can
    specify the name of the quantity as an item::

         >>> grid['density']

    which is no longer a :class:`~hyperion.grid.VoronoiGrid` object, but
    a :class:`~hyperion.grid.VoronoiGridView` object. When setting
    this for the first time, this can be set either to another
    :class:`~hyperion.grid.VoronoiGridView` object, an external h5py
    link, or an empty list. For example, the following should work:

        >>> grid['density_new'] = grid['density']

    :class:`~hyperion.grid.VoronoiGridView` objects allow the
    specific dust population to be selected as an index:

        >>> grid['density'][0]

    Which is also a :class:`~hyperion.grid.VoronoiGridView` object. The
    data can then be accessed with the ``array`` attribute::

        >>> grid['density'][0].array

    which is a 1-d array of the requested quantity.
    '''

    def __init__(self, *args, **kwargs):

        self.shape = None

        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None

        self._x = None
        self._y = None
        self._z = None

        self.quantities = {}

        self.voronoi_table = None

        try:
            self._verbose = kwargs.pop('verbose')
        except KeyError:
            self._verbose = False

        self._freeze()

        if len(args) > 0:
            if isinstance(args[0], VoronoiGrid):
                self.set_points(args[0].x, args[0].y, args[0].z,
                                xmin=args[0].xmin, xmax=args[0].xmax,
                                ymin=args[0].ymin, ymax=args[0].ymax,
                                zmin=args[0].zmin, zmax=args[0].zmax)
                self.voronoi_table = args[0].voronoi_table
            else:
                self.set_points(*args, **kwargs)

    def set_points(self, x, y, z,
                   xmin=None, xmax=None,
                   ymin=None, ymax=None,
                   zmin=None, zmax=None):

        if type(x) in [list, tuple]:
            x = np.array(x)
        if type(y) in [list, tuple]:
            y = np.array(y)
        if type(z) in [list, tuple]:
            z = np.array(z)

        if not is_numpy_array(x) or x.ndim != 1:
            raise ValueError("x should be a 1-D sequence")
        if not is_numpy_array(y) or y.ndim != 1:
            raise ValueError("y should be a 1-D sequence")
        if not is_numpy_array(z) or z.ndim != 1:
            raise ValueError("z should be a 1-D sequence")

        # Make sure we cast to 64-bit floating point values otherwise strange
        # things will happen in the fortran code.
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        # Find grid shape
        self.shape = (len(x),)

        # If limits were specified, use those

        bounds = [xmin, xmax, ymin, ymax, zmin, zmax]

        if None in bounds:

            if bounds.count(None) != 6:
                raise ValueError("Either all or no limits should be specified")

            # Add 1% border around points
            delta_x = x.max() - x.min()
            delta_y = y.max() - y.min()
            delta_z = z.max() - z.min()
            self.xmin = x.min() - delta_x / 100.
            self.xmax = x.max() + delta_x / 100.
            self.ymin = y.min() - delta_y / 100.
            self.ymax = y.max() + delta_y / 100.
            self.zmin = z.min() - delta_z / 100.
            self.zmax = z.max() + delta_z / 100.

        else:

            self.xmin = float(xmin)
            self.xmax = float(xmax)
            self.ymin = float(ymin)
            self.ymax = float(ymax)
            self.zmin = float(zmin)
            self.zmax = float(zmax)

        self._x = x
        self._y = y
        self._z = z

    @property
    def voronoi_table(self):

        if self._voronoi_table is None or self._voronoi_table.meta['geometry'].decode('utf-8') != self.get_geometry_id():

            from .voronoi_helpers import voronoi_grid

            logger.info("Updating Voronoi Tesselation")

            # Compute the Voronoi tesselation
            points = np.array([self._x, self._y, self._z]).transpose()
            mesh = voronoi_grid(points, np.array([[self.xmin, self.xmax],
                                                  [self.ymin, self.ymax],
                                                  [self.zmin, self.zmax]]),
                                verbose=self._verbose)
            self._voronoi_table = mesh.neighbours_table
            self._voronoi_table.meta['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        return self._voronoi_table

    @voronoi_table.setter
    def voronoi_table(self, table):
        if table is None or table.meta['geometry'].decode('utf-8') == self.get_geometry_id():
            self._voronoi_table = table
        else:
            raise ValueError("geometry does not match: expected {0} but got {1}".format(self.get_geometry_id(), table.meta['geometry']))

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    def __getattr__(self, attribute):
        if attribute == 'n_dust':
            n_dust = None
            for quantity in self.quantities:
                n_dust_q, shape_q = single_grid_dims(self.quantities[quantity], ndim=1)
                if n_dust is None:
                    n_dust = n_dust_q
                elif n_dust_q is not None:
                    if n_dust != n_dust_q:
                        raise ValueError("Not all dust lists in the grid have the same size")
            return n_dust
        else:
            return FreezableClass.__getattribute__(self, attribute)

    def _check_array_dimensions(self, array=None):
        '''
        Check that a grid's array dimensions agree with this grid's metadata

        Parameters
        ----------
        array : np.ndarray or list of np.ndarray, optional
            The array for which to test the dimensions. If this is not
            specified, this method performs a self-consistency check of array
            dimensions and meta-data.
        '''

        n_pop_ref = None

        if isinstance(array, VoronoiGridView):
            array = array.quantities[array.viewed_quantity]

        for quantity in self.quantities:

            if array is None:
                n_pop, shape = single_grid_dims(self.quantities[quantity], ndim=1)
            else:
                n_pop, shape = single_grid_dims(array, ndim=1)

            if shape != self.shape:
                raise ValueError("Quantity arrays do not have the right "
                                 "dimensions: %s instead of %s"
                                 % (shape, self.shape))

            if n_pop is not None:
                if n_pop_ref is None:
                    n_pop_ref = n_pop
                elif n_pop != n_pop_ref:
                    raise ValueError("Not all dust lists in the grid have the same size")

    def read(self, group, quantities='all'):
        '''
        Read the geometry and physical quantities from an voronoi grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid from. This group should contain
            groups named 'Geometry' and 'Quantities'.
        quantities : 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Read in geometry
        self.read_geometry(group['Geometry'])

        # Read in physical quantities
        self.read_quantities(group['Quantities'], quantities=quantities)

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

    def read_geometry(self, group):
        '''
        Read in geometry information from a cartesian grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid geometry from.
        '''

        if group.attrs['grid_type'].decode('utf-8') != 'vor':
            raise ValueError("Grid is not an voronoi")

        coords = group['cells']['coordinates']

        self.set_points(coords[:,0], coords[:,1], coords[:,2],
                        xmin=group.attrs['xmin'], xmax=group.attrs['xmax'],
                        ymin=group.attrs['ymin'], ymax=group.attrs['ymax'],
                        zmin=group.attrs['zmin'], zmax=group.attrs['zmax'])

        # Check that advertised hash matches real hash
        if group.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

        # Avoid re-computing Voronoi table
        self.voronoi_table = Table.read(group['cells'], format='hdf5')

    def read_quantities(self, group, quantities='all'):
        '''
        Read in physical quantities from a cartesian grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid quantities from
        quantities : 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Read in physical quantities
        if quantities is not None:
            for quantity in group:
                if quantities == 'all' or quantity in quantities:
                    array = np.array(group[quantity])
                    if array.ndim == 2:  # if array is 2D, it is a list of 1D arrays
                        self.quantities[quantity] = [array[i] for i in range(array.shape[0])]
                    else:
                        self.quantities[quantity] = array

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

    @property
    def volumes(self):
        return np.array(self.voronoi_table['volume'])

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the voronoi grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to write the grid to
        quantities : 'all' or list
            Which physical quantities to write out. Use 'all' to write out all
            quantities or a list of strings to write only specific quantities.
        copy : bool
            Whether to copy external links, or leave them as links.
        absolute_paths : bool
            If copy is False, then this indicates whether to use absolute or
            relative paths for links.
        compression : bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype : type
            The datatype to use to write the wall positions (ignored for this kind of grid)
        physics_dtype : type
            The datatype to use to write the physical quantities
        '''

        # Create HDF5 groups if needed

        if 'Geometry' not in group:
            g_geometry = group.create_group('Geometry')
        else:
            g_geometry = group['Geometry']

        if 'Quantities' not in group:
            g_quantities = group.create_group('Quantities')
        else:
            g_quantities = group['Quantities']

        # Write out geometry

        g_geometry.attrs['grid_type'] = np.string_('vor'.encode('utf-8'))
        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        g_geometry.attrs['xmin'] = self.xmin
        g_geometry.attrs['xmax'] = self.xmax
        g_geometry.attrs['ymin'] = self.ymin
        g_geometry.attrs['ymax'] = self.ymax
        g_geometry.attrs['zmin'] = self.zmin
        g_geometry.attrs['zmax'] = self.zmax

        # Compute Voronoi tesselation
        voronoi_table = self.voronoi_table

        # Fix nan/zero
        voronoi_table['volume'][voronoi_table['volume'] <= 0.] = -1.
        voronoi_table['volume'][np.isnan(voronoi_table['volume'])] = -1.
        voronoi_table['volume'][np.isinf(voronoi_table['volume'])] = -1.

        if voronoi_table['neighbours'].dtype != np.int32:
            raise TypeError("neighbours should be int32")

        voronoi_table.write(g_geometry, path="cells", compression=True)

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

        # Write out physical quantities

        for quantity in self.quantities:
            if quantities == 'all' or quantity in quantities:
                if isinstance(self.quantities[quantity], h5py.ExternalLink):
                    link_or_copy(g_quantities, quantity, self.quantities[quantity], copy, absolute_paths=absolute_paths)
                else:
                    dset = g_quantities.create_dataset(quantity, data=self.quantities[quantity],
                                                       compression=compression,
                                                       dtype=physics_dtype)
                    dset.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

    def write_single_array(self, group, name, array, copy=True, absolute_paths=False, compression=True, physics_dtype=float):
        '''
        Write out a single quantity, checking for consistency with geometry

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to write the grid to
        name : str
            The name of the array in the group
        array : np.ndarray
            The array to write out
        copy : bool
            Whether to copy external links, or leave them as links.
        absolute_paths : bool
            If copy is False, then this indicates whether to use absolute or
            relative paths for links.
        compression : bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype : type
            The datatype to use to write the wall positions
        physics_dtype : type
            The datatype to use to write the physical quantities
        '''

        # Check consistency of array dimensions with grid
        self._check_array_dimensions(array)

        if isinstance(array, h5py.ExternalLink):
            link_or_copy(group, name, array, copy, absolute_paths=absolute_paths)
        else:
            dset = group.create_dataset(name, data=array,
                                        compression=compression,
                                        dtype=physics_dtype)
            dset.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

    def get_geometry_id(self):
        # The grid is uniquely defined by the points and the bounds
        geo_hash = hashlib.md5()
        geo_hash.update(self.x.tostring())
        geo_hash.update(self.y.tostring())
        geo_hash.update(self.z.tostring())
        geo_hash.update(struct.pack('>d', self.xmin))
        geo_hash.update(struct.pack('>d', self.xmax))
        geo_hash.update(struct.pack('>d', self.ymin))
        geo_hash.update(struct.pack('>d', self.ymax))
        geo_hash.update(struct.pack('>d', self.zmin))
        geo_hash.update(struct.pack('>d', self.zmax))
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return VoronoiGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, VoronoiGridView):
            if self.refined is None:
                logger.warn("No geometry in target grid - copying from original grid")
                self.set_points(value.x, value.y, value.z)
            self.quantities[item] = deepcopy(value.quantities[value.viewed_quantity])
        elif isinstance(value, h5py.ExternalLink):
            self.quantities[item] = value
        elif value == []:
            self.quantities[item] = []
        else:
            raise ValueError('value should be an empty list, and ExternalLink, or a VoronoiGridView instance')

    def __contains__(self, item):
        return self.quantities.__contains__(item)

    def reset_quantities(self):
        self.quantities = {}

    def add_derived_quantity(self, name, function):
        if name in self.quantities:
            raise KeyError(name + ' already exists')
        function(self.quantities)


class VoronoiGridView(VoronoiGrid):

    def __init__(self, grid, quantity):
        self.viewed_quantity = quantity
        VoronoiGrid.__init__(self)
        self.set_points(grid.x, grid.y, grid.z,
                        xmin=grid.xmin, xmax=grid.xmax,
                        ymin=grid.ymin, ymax=grid.ymax,
                        zmin=grid.zmin, zmax=grid.zmax)
        self.voronoi_table = grid.voronoi_table
        self.quantities = {quantity: grid.quantities[quantity]}

    def append(self, grid):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        grid : 1D Numpy array or VoronoiGridView instance
            The grid to copy the quantity from
        '''
        if isinstance(grid, VoronoiGridView):
            if self.quantities[self.viewed_quantity] is grid.quantities[grid.viewed_quantity]:
                raise Exception("Calling append recursively")
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("Can only append a single grid")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity].append(deepcopy(grid.quantities[grid.viewed_quantity]))
        elif isinstance(grid, np.ndarray):
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity].append(deepcopy(grid))
        else:
            raise ValueError("grid should be a Numpy array or an VoronoiGridView instance")

    def add(self, grid):
        '''
        Used to add quantities from another grid

        Parameters
        ----------
        grid : 1D Numpy array or VoronoiGridView instance
            The grid to copy the quantity from
        '''
        if type(self.quantities[self.viewed_quantity]) is list:
            raise Exception("need to first specify the item to add to")
        if isinstance(grid, VoronoiGridView):
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("need to first specify the item to add")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity] += grid.quantities[grid.viewed_quantity]
        elif isinstance(grid, np.ndarray):
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity] += grid
        else:
            raise ValueError("grid should be a Numpy array or an VoronoiGridView instance")

    def __getitem__(self, item):
        if type(item) is int:
            grid = VoronoiGridView(self, self.viewed_quantity)
            grid.quantities = {grid.viewed_quantity: grid.quantities[grid.viewed_quantity][item]}
            return grid
        else:
            return VoronoiGrid.__getitem__(self, item)

    def __getattr__(self, attribute):
        if attribute == 'array':
            return self.quantities[self.viewed_quantity]
        else:
            return VoronoiGrid.__getattr__(self, attribute)
