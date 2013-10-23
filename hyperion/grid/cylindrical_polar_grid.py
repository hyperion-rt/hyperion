from __future__ import print_function, division

import hashlib
from copy import deepcopy

import h5py
import numpy as np

from ..util.meshgrid import meshgrid_nd
from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing, link_or_copy
from astropy import log as logger
from .grid_helpers import single_grid_dims


class CylindricalPolarGrid(FreezableClass):
    '''
    A cylindrical polar grid.

    The grid can be initialized by passing the w, z, and phi coordinates of cell walls::

        >>> grid = CylindricalPolarGrid(w_wall, z_wall, p_wall)

    where ``w_wall``, ``z_wall``, and ``p_wall`` are 1-d sequences of wall
    positions. The number of cells in the resulting grid will be one less
    in each dimension that the length of these arrays.

    :class:`~hyperion.grid.CylindricalPolarGrid` objects may contain multiple
    quantities (e.g. density, specific energy). To access these, you can
    specify the name of the quantity as an item::

         >>> grid['density']

    which is no longer a :class:`~hyperion.grid.CylindricalPolarGrid` object, but
    a :class:`~hyperion.grid.CylindricalPolarGridView` object. When setting
    this for the first time, this can be set either to another
    :class:`~hyperion.grid.CylindricalPolarGridView` object, an external h5py
    link, or an empty list. For example, the following should work:

        >>> grid['density_new'] = grid['density']

    :class:`~hyperion.grid.CylindricalPolarGridView` objects allow the
    specific dust population to be selected as an index:

        >>> grid['density'][0]

    Which is also a :class:`~hyperion.grid.CylindricalPolarGridView` object. The
    data can then be accessed with the ``array`` attribute::

        >>> grid['density'][0].array

    which is a 3-d array of the requested quantity.
    '''

    def __init__(self, *args):

        self.shape = None

        self.w_wall = None
        self.z_wall = None
        self.p_wall = None

        self.w = None
        self.z = None
        self.p = None

        self.gw = None
        self.gz = None
        self.gp = None

        self.volumes = None
        self.areas = None
        self.widths = None

        self.quantities = {}

        self._freeze()

        if len(args) > 0:
            if isinstance(args[0], CylindricalPolarGrid):
                self.set_walls(args[0].w_wall, args[0].z_wall, args[0].p_wall)
            else:
                self.set_walls(*args)

    def set_walls(self, w_wall, z_wall, p_wall):

        if type(w_wall) in [list, tuple]:
            w_wall = np.array(w_wall)
        if type(z_wall) in [list, tuple]:
            z_wall = np.array(z_wall)
        if type(p_wall) in [list, tuple]:
            p_wall = np.array(p_wall)

        if not is_numpy_array(w_wall) or w_wall.ndim != 1:
            raise ValueError("w_wall should be a 1-D sequence")
        if not is_numpy_array(z_wall) or z_wall.ndim != 1:
            raise ValueError("z_wall should be a 1-D sequence")
        if not is_numpy_array(p_wall) or p_wall.ndim != 1:
            raise ValueError("p_wall should be a 1-D sequence")

        if not monotonically_increasing(w_wall):
            raise ValueError("w_wall should be monotonically increasing")
        if not monotonically_increasing(z_wall):
            raise ValueError("z_wall should be monotonically increasing")
        if not monotonically_increasing(p_wall):
            raise ValueError("p_wall should be monotonically increasing")

        if np.any(p_wall < 0.) or np.any(p_wall > 2. * np.pi):
            raise ValueError("p_wall values be in the range [0:2*pi]")

        # Find grid shape
        self.shape = (len(p_wall) - 1, len(z_wall) - 1, len(w_wall) - 1)

        # Store wall positions
        self.w_wall = w_wall
        self.z_wall = z_wall
        self.p_wall = p_wall

        # Compute cell centers
        if w_wall[0] == 0.:
            self.w = np.zeros(len(w_wall) - 1)
            self.w[0] = w_wall[1] / 2.
            self.w[1:] = 10. ** ((np.log10(w_wall[1:-1]) + np.log10(w_wall[2:])) / 2.)
        else:
            self.w = 10. ** ((np.log10(w_wall[:-1]) + np.log10(w_wall[1:])) / 2.)

        self.z = (z_wall[:-1] + z_wall[1:]) / 2.
        self.p = (p_wall[:-1] + p_wall[1:]) / 2.

        # Generate 3D versions of r, t, p
        #(each array is 3D and defined in every cell)
        self.gw, self.gz, self.gp = meshgrid_nd(self.w, self.z, self.p)

        # Generate 3D versions of the inner and outer wall positions respectively
        gw_wall_min, gz_wall_min, gp_wall_min = \
                    meshgrid_nd(w_wall[:-1], z_wall[:-1], p_wall[:-1])

        gw_wall_max, gz_wall_max, gp_wall_max = \
                    meshgrid_nd(w_wall[1:], z_wall[1:], p_wall[1:])

        # USEFUL QUANTITIES

        dr = gw_wall_max - gw_wall_min
        dr2 = gw_wall_max ** 2 - gw_wall_min ** 2
        dz = gz_wall_max - gz_wall_min
        dp = gp_wall_max - gp_wall_min

        # CELL VOLUMES

        #   dV = dr * dz * (r*dphi)
        #    V = [r_2^2 - r_1^2] / 2. * [z_2 - z_1] * [phi_2 - phi_1]

        self.volumes = dr2 * dz * dp / 2.

        # WALL AREAS

        self.areas = np.zeros((6,) + self.shape)

        # R walls:
        #   dA = r * dz * dphi
        #    A = r * [z 2 - z_1] * [phi_2 - phi_1]

        self.areas[0, :, :, :] = gw_wall_min * dz * dp
        self.areas[1, :, :, :] = gw_wall_max * dz * dp

        # z walls:
        #   dA = r * dr * dphi
        #    A = 0.5 * [r_2^2 - r_1^2] * [phi_2 - phi_1]

        self.areas[2, :, :, :] = 0.5 * dr2 * dp
        self.areas[3, :, :, :] = 0.5 * dr2 * dp

        # Phi walls:
        #   dA = dr * dz
        #    A = [r_2 - r_1] * [z_2 - z_1]

        self.areas[4, :, :, :] = dr * dz
        self.areas[5, :, :, :] = dr * dz

        # CELL WIDTHS

        self.widths = np.zeros((3,) + self.shape)

        # R direction:
        #   dS = dr
        #    S = r_2 - r_1

        self.widths[0, :, :, :] = dr

        # z direction:
        #   dS = dz
        #    S = [z_2 - z_1]

        self.widths[1, :, :, :] = dz

        # Phi direction:
        #   dS = r * dphi
        #    S = r * [phi_2 - phi_1]

        self.widths[2, :, :, :] = self.gw * dp

    def __getattr__(self, attribute):
        if attribute == 'n_dust':
            n_dust = None
            for quantity in self.quantities:
                n_dust_q, shape_q = single_grid_dims(self.quantities[quantity])
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

        if isinstance(array, CylindricalPolarGridView):
            array = array.quantities[array.viewed_quantity]

        for quantity in self.quantities:

            if array is None:
                n_pop, shape = single_grid_dims(self.quantities[quantity])
            else:
                n_pop, shape = single_grid_dims(array)

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
        Read the geometry and physical quantities from a cylindrical polar grid

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
        Read in geometry information from a cylindrical polar grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid geometry from.
        '''

        if group.attrs['grid_type'].decode('utf-8') != 'cyl_pol':
            raise ValueError("Grid is not cylindrical polar")

        self.set_walls(group['walls_1']['w'],
                       group['walls_2']['z'],
                       group['walls_3']['p'])

        # Check that advertised hash matches real hash
        if group.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

    def read_quantities(self, group, quantities='all'):
        '''
        Read in physical quantities from a cylindrical polar grid

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
                    if array.ndim == 4:  # if array is 4D, it is a list of 3D arrays
                        self.quantities[quantity] = [array[i] for i in range(array.shape[0])]
                    else:
                        self.quantities[quantity] = array

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the cylindrical polar grid

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
            The datatype to use to write the wall positions
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

        g_geometry.attrs['grid_type'] = np.string_('cyl_pol'.encode('utf-8'))
        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        dset = g_geometry.create_dataset("walls_1", data=np.array(list(zip(self.w_wall)), dtype=[('w', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = np.string_('cm'.encode('utf-8'))

        dset = g_geometry.create_dataset("walls_2", data=np.array(list(zip(self.z_wall)), dtype=[('z', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = np.string_('cm'.encode('utf-8'))

        dset = g_geometry.create_dataset("walls_3", data=np.array(list(zip(self.p_wall)), dtype=[('p', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = np.string_('rad'.encode('utf-8'))

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
        geo_hash.update(self.w_wall.tostring())
        geo_hash.update(self.z_wall.tostring())
        geo_hash.update(self.p_wall.tostring())
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return CylindricalPolarGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, CylindricalPolarGridView):
            if self.w_wall is None and self.z_wall is None and self.p_wall is None:
                logger.warn("No geometry in target grid - copying from original grid")
                self.set_walls(value.w_wall, value.z_wall, value.p_wall)
            self.quantities[item] = deepcopy(value.quantities[value.viewed_quantity])
        elif isinstance(value, h5py.ExternalLink):
            self.quantities[item] = value
        elif value == []:
            self.quantities[item] = []
        else:
            raise ValueError('value should be an empty list, and ExternalLink, or a CylindricalPolarGridView instance')

    def __contains__(self, item):
        return self.quantities.__contains__(item)

    def reset_quantities(self):
        self.quantities = {}

    def add_derived_quantity(self, name, function):
        if name in self.quantities:
            raise KeyError(name + ' already exists')
        function(self.quantities)


class CylindricalPolarGridView(CylindricalPolarGrid):

    def __init__(self, grid, quantity):
        self.viewed_quantity = quantity
        CylindricalPolarGrid.__init__(self)
        self.set_walls(grid.w_wall, grid.z_wall, grid.p_wall)
        self.quantities = {quantity: grid.quantities[quantity]}

    def append(self, grid):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        grid : 3D Numpy array or CylindricalPolarGridView instance
            The grid to copy the quantity from
        '''
        if isinstance(grid, CylindricalPolarGridView):
            if self.quantities[self.viewed_quantity] is grid.quantities[grid.viewed_quantity]:
                raise Exception("Calling append recursively")
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("Can only append a single grid")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity].append(deepcopy(grid.quantities[grid.viewed_quantity]))
        elif type(grid) is np.ndarray:
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity].append(deepcopy(grid))
        else:
            raise ValueError("grid should be a Numpy array or a CylindricalPolarGridView instance")

    def add(self, grid):
        '''
        Used to add quantities from another grid

        Parameters
        ----------
        grid : 3D Numpy array or CylindricalPolarGridView instance
            The grid to copy the quantity from
        '''
        if type(self.quantities[self.viewed_quantity]) is list:
            raise Exception("need to first specify the item to add to")
        if isinstance(grid, CylindricalPolarGridView):
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("need to first specify the item to add")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity] += grid.quantities[grid.viewed_quantity]
        elif type(grid) is np.ndarray:
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity] += grid
        else:
            raise ValueError("grid should be a Numpy array or a CylindricalPolarGridView instance")

    def __getitem__(self, item):
        if type(item) is int:
            grid = CylindricalPolarGridView(self, self.viewed_quantity)
            grid.quantities = {grid.viewed_quantity: grid.quantities[grid.viewed_quantity][item]}
            return grid
        else:
            return CylindricalPolarGrid.__getitem__(self, item)

    def __getattr__(self, attribute):
        if attribute == 'array':
            return self.quantities[self.viewed_quantity]
        else:
            return CylindricalPolarGrid.__getattr__(self, attribute)
