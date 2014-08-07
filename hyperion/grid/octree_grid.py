from __future__ import print_function, division

import struct
import hashlib
from copy import deepcopy

import h5py
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, link_or_copy
from astropy import log as logger
from .grid_helpers import single_grid_dims


class OctreeGrid(FreezableClass):
    '''
    An octree grid.

    To initialize an Octree object, use::

        >>> grid = OctreeGrid(x, y, z, dx, dy, dz, refined)

    where ``x``, ``y``, and ``z`` are the cartesian coordinates of the center
    of the grid, ``dx``, ``dy``, and ``dz`` are the half-widths of the grid,
    and ``refined`` is a sequence of boolean values that indicate whether a
    given cell is refined.

    The first value of the ``refined`` sequence indicates whether the parent
    cell is sub-divided. If it is, then the the second element indicates
    whether the first cell of the parent cell is sub-divided. If it isn't,
    then the next value indicates whether the second cell of the parent cell
    is sub-divided. If it is, then we need to specify the booleans for all the
    children of that cell before we move to the third cell of the parent cell.

    For example, the simplest grid is a single cell that is not sub-divided::

        refined = [False]

    The next simplest grid is a single grid cell that is only sub-divided once::

        refined = [True, False, False, False, False, False, False, False, False]

    It is easier to picture this as a hierarchy::

        refined = [True,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     ]

    If we sub-divide the third sub-cell in the parent cell into cells that are
    themselves not sub-divided, we get::

        refined = [True,
                     False,
                     False,
                     True,
                       False,
                       False,
                       False,
                       False,
                       False,
                       False,
                       False,
                       False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     ]

    and so on. The order of the sub-cells is first along x, then along y, then
    along z.

    :class:`~hyperion.grid.OctreeGrid` objects may contain multiple
    quantities (e.g. density, specific energy). To access these, you can
    specify the name of the quantity as an item::

         >>> grid['density']

    which is no longer an :class:`~hyperion.grid.OctreeGrid` object, but
    a :class:`~hyperion.grid.OctreeGridView` object. When setting
    this for the first time, this can be set either to another
    :class:`~hyperion.grid.OctreeGridView` object, an external h5py
    link, or an empty list. For example, the following should work:

        >>> grid['density_new'] = grid['density']

    :class:`~hyperion.grid.OctreeGridView` objects allow the
    specific dust population to be selected as an index:

        >>> grid['density'][0]

    Which is also an :class:`~hyperion.grid.OctreeGridView` object.
    '''

    _validate_cache = {}

    def __init__(self, *args):

        self.shape = None

        self.x = None
        self.y = None
        self.z = None

        self.dx = None
        self.dy = None
        self.dz = None

        self.refined = None

        self.quantities = {}

        self._freeze()

        if len(args) > 0:
            if isinstance(args[0], OctreeGrid):
                self.set_walls(args[0].x, args[0].y, args[0].z,
                               args[0].dx, args[0].dy, args[0].dz,
                               args[0].refined)
            else:
                self.set_walls(*args)

    def set_walls(self, x, y, z, dx, dy, dz, refined):

        if not np.isscalar(x):
            raise ValueError("x should be a scalar value")
        if not np.isreal(x):
            raise ValueError("x should be a numerical value")

        if not np.isscalar(y):
            raise ValueError("y should be a scalar value")
        if not np.isreal(y):
            raise ValueError("y should be a numerical value")

        if not np.isscalar(z):
            raise ValueError("z should be a scalar value")
        if not np.isreal(z):
            raise ValueError("z should be a numerical value")

        if not np.isscalar(dx):
            raise ValueError("dx should be a scalar value")
        if not np.isreal(dx):
            raise ValueError("dx should be a numerical value")

        if not np.isscalar(dy):
            raise ValueError("dy should be a scalar value")
        if not np.isreal(dy):
            raise ValueError("dy should be a numerical value")

        if not np.isscalar(dz):
            raise ValueError("dz should be a scalar value")
        if not np.isreal(dz):
            raise ValueError("dz should be a numerical value")

        if type(refined) in [list, tuple]:
            refined = np.array(refined)

        if not is_numpy_array(refined) or refined.ndim != 1 or refined.dtype != bool:
            raise ValueError("refined should be a 1-D boolean sequence")

        self.x = x
        self.y = y
        self.z = z

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.refined = refined

        self.shape = (len(refined),)

    @property
    def refined(self):
        return self._refined

    @refined.setter
    def refined(self, value):

        if value is None:
            self._refined = None
            return

        if not (len(value) - 1) % 8 == 0:
            raise ValueError("refined should have shape 8 * n + 1")

        self._refined = self._validate(value)

    def _validate(self, value):

        value_hash = hashlib.md5(value.tostring()).hexdigest()

        if value_hash in self._validate_cache:
            return value

        logger.info("Checking consistency of refined array")

        # Check that refined array reduces to a single False if removing all
        # levels of refinement.
        refined_str = value.tostring()
        previous = ''
        while True:
            refined_str = refined_str.replace(b'\x01\x00\x00\x00\x00\x00\x00\x00\x00', b'\x00')
            if refined_str == previous:
                break
            else:
                previous = refined_str
        if refined_str != b'\x00':
            raise ValueError('refined should reduce to a single False value if removing levels in hierarchy')

        def check_recursive(refined, current_i=0, max_level=0):
            if refined[current_i]:
                current_i += 1
                max_levels = []
                for i in range(8):
                    current_i, max_level_indiv = check_recursive(refined, current_i, max_level+1)
                    max_levels.append(max_level_indiv)
                max_level = max(max_levels)
            else:
                current_i += 1
            return current_i, max_level

        try:
            final_i, max_level = check_recursive(value)
        except IndexError:
            raise ValueError("refined array is not self-consistent")

        logger.info("Setting refined with maximum depth of {0} levels".format(max_level))

        if max_level > 20:
            logger.warn("Number of levels in octree is high ({0})".format(max_level))

        self._validate_cache[value_hash] = True

        return value

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

    @property
    def limits(self):
        from hyperion.importers._discretize_sph import _get_positions_widths
        xc, yc, zc, xw, yw, zw = _get_positions_widths(self.refined,
                                                       self.x, self.y, self.z,
                                                       self.dx, self.dy, self.dz)
        return xc - xw, xc + xw, yc - yw, yc + yw, zc - zw, zc + zw

    @property
    def volumes(self):
        """
        The volumes of all the cells in the octree
        """

        xmin, xmax, ymin, ymax, zmin, zmax = self.limits
        return (xmax - xmin) * (ymax - ymin) * (zmax - zmin)

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

        if isinstance(array, OctreeGridView):
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
        Read the geometry and physical quantities from an octree grid

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

        if group.attrs['grid_type'].decode('utf-8') != 'oct':
            raise ValueError("Grid is not an octree")

        self.set_walls(group.attrs['x'],
                       group.attrs['y'],
                       group.attrs['z'],
                       group.attrs['dx'],
                       group.attrs['dy'],
                       group.attrs['dz'],
                       group['cells']['refined'].astype(bool))

        # Check that advertised hash matches real hash
        if group.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

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

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the octree grid

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

        g_geometry.attrs['grid_type'] = np.string_('oct'.encode('utf-8'))
        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        g_geometry.attrs['x'] = self.x
        g_geometry.attrs['y'] = self.y
        g_geometry.attrs['z'] = self.z
        g_geometry.attrs['dx'] = self.dx
        g_geometry.attrs['dy'] = self.dy
        g_geometry.attrs['dz'] = self.dz

        dset = g_geometry.create_dataset("cells", data=np.array(list(zip(self.refined)), dtype=[('refined', np.int32)]), compression=compression)

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
        geo_hash = hashlib.md5()
        geo_hash.update(struct.pack('>d', self.x))
        geo_hash.update(struct.pack('>d', self.y))
        geo_hash.update(struct.pack('>d', self.z))
        geo_hash.update(struct.pack('>d', self.dx))
        geo_hash.update(struct.pack('>d', self.dy))
        geo_hash.update(struct.pack('>d', self.dz))
        geo_hash.update(self.refined.tostring())
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return OctreeGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, OctreeGridView):
            if self.refined is None:
                logger.warn("No geometry in target grid - copying from original grid")
                self.set_walls(value.x, value.y, value.z, value.dx, value.dy, value.dz, value.refined)
            self.quantities[item] = deepcopy(value.quantities[value.viewed_quantity])
        elif isinstance(value, h5py.ExternalLink):
            self.quantities[item] = value
        elif value == []:
            self.quantities[item] = []
        else:
            raise ValueError('value should be an empty list, and ExternalLink, or a OctreeGridView instance')

    def __contains__(self, item):
        return self.quantities.__contains__(item)

    def reset_quantities(self):
        self.quantities = {}

    def add_derived_quantity(self, name, function):
        if name in self.quantities:
            raise KeyError(name + ' already exists')
        function(self.quantities)

    def to_yt(self, dust_id=0):
        '''
        Convert AMR grid to a yt object (requires yt)

        Parameters
        ----------
        dust_id : int, optional
            The ID of the dust population to extract. If not set, this
            defaults to 0 (the first dust population).
        '''
        from yt_wrappers import octree_grid_to_yt_stream
        return octree_grid_to_yt_stream(self, dust_id)


class OctreeGridView(OctreeGrid):

    def __init__(self, grid, quantity):
        self.viewed_quantity = quantity
        OctreeGrid.__init__(self)
        self.set_walls(grid.x, grid.y, grid.z, grid.dx, grid.dy, grid.dz, grid.refined)
        self.quantities = {quantity: grid.quantities[quantity]}

    def append(self, grid):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        grid : 1D Numpy array or OctreeGridView instance
            The grid to copy the quantity from
        '''
        if isinstance(grid, OctreeGridView):
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
            raise ValueError("grid should be a Numpy array or an OctreeGridView instance")

    def add(self, grid):
        '''
        Used to add quantities from another grid

        Parameters
        ----------
        grid : 1D Numpy array or OctreeGridView instance
            The grid to copy the quantity from
        '''
        if type(self.quantities[self.viewed_quantity]) is list:
            raise Exception("need to first specify the item to add to")
        if isinstance(grid, OctreeGridView):
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("need to first specify the item to add")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity] += grid.quantities[grid.viewed_quantity]
        elif isinstance(grid, np.ndarray):
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity] += grid
        else:
            raise ValueError("grid should be a Numpy array or an OctreeGridView instance")

    def __getitem__(self, item):
        if type(item) is int:
            grid = OctreeGridView(self, self.viewed_quantity)
            grid.quantities = {grid.viewed_quantity: grid.quantities[grid.viewed_quantity][item]}
            return grid
        else:
            return OctreeGrid.__getitem__(self, item)

    def __getattr__(self, attribute):
        if attribute == 'array':
            return self.quantities[self.viewed_quantity]
        else:
            return OctreeGrid.__getattr__(self, attribute)
