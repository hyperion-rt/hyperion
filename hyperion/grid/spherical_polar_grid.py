from __future__ import print_function, division

import hashlib
from copy import deepcopy

import h5py
import numpy as np

from ..util.meshgrid import meshgrid_nd
from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing, link_or_copy
from astropy import log as logger
from .grid_helpers import single_grid_dims


class SphericalPolarGrid(FreezableClass):
    '''
    A spherical polar grid.

    The grid can be initialized by passing the r, theta, and phi coordinates of cell walls::

        >>> grid = SphericalPolarGrid(r_wall, t_wall, p_wall)

    where ``r_wall``, ``t_wall``, and ``p_wall`` are 1-d sequences of wall
    positions. The number of cells in the resulting grid will be one less
    in each dimension that the length of these arrays.

    :class:`~hyperion.grid.SphericalPolarGrid` objects may contain multiple
    quantities (e.g. density, specific energy). To access these, you can
    specify the name of the quantity as an item::

         >>> grid['density']

    which is no longer a :class:`~hyperion.grid.SphericalPolarGrid` object, but
    a :class:`~hyperion.grid.SphericalPolarGridView` object. When setting
    this for the first time, this can be set either to another
    :class:`~hyperion.grid.SphericalPolarGridView` object, an external h5py
    link, or an empty list. For example, the following should work:

        >>> grid['density_new'] = grid['density']

    :class:`~hyperion.grid.SphericalPolarGridView` objects allow the
    specific dust population to be selected as an index:

        >>> grid['density'][0]

    Which is also a :class:`~hyperion.grid.SphericalPolarGridView` object. The
    data can then be accessed with the ``array`` attribute::

        >>> grid['density'][0].array

    which is a 3-d array of the requested quantity.
    '''

    def __init__(self, *args):

        self.shape = None

        self.r_wall = None
        self.t_wall = None
        self.p_wall = None

        self.r = None
        self.t = None
        self.p = None

        self.gr = None
        self.gt = None
        self.gp = None

        self.gw = None
        self.gz = None

        self.volumes = None
        self.areas = None
        self.widths = None

        self.quantities = {}

        self._freeze()

        if len(args) > 0:
            if isinstance(args[0], SphericalPolarGrid):
                self.set_walls(args[0].r_wall, args[0].t_wall, args[0].p_wall)
            else:
                self.set_walls(*args)

    def set_walls(self, r_wall, t_wall, p_wall):

        if type(r_wall) in [list, tuple]:
            r_wall = np.array(r_wall)
        if type(t_wall) in [list, tuple]:
            t_wall = np.array(t_wall)
        if type(p_wall) in [list, tuple]:
            p_wall = np.array(p_wall)

        if not is_numpy_array(r_wall) or r_wall.ndim != 1:
            raise ValueError("r_wall should be a 1-D sequence")
        if not is_numpy_array(t_wall) or t_wall.ndim != 1:
            raise ValueError("t_wall should be a 1-D sequence")
        if not is_numpy_array(p_wall) or p_wall.ndim != 1:
            raise ValueError("p_wall should be a 1-D sequence")

        if not monotonically_increasing(r_wall):
            raise ValueError("r_wall should be monotonically increasing")
        if not monotonically_increasing(t_wall):
            raise ValueError("t_wall should be monotonically increasing")
        if not monotonically_increasing(p_wall):
            raise ValueError("p_wall should be monotonically increasing")

        if np.any(t_wall < 0.) or np.any(t_wall > np.pi):
            raise ValueError("t_wall values be in the range [0:pi]")
        if np.any(p_wall < 0.) or np.any(p_wall > 2. * np.pi):
            raise ValueError("p_wall values be in the range [0:2*pi]")

        # Find number of grid cells
        self.shape = (len(p_wall) - 1, len(t_wall) - 1, len(r_wall) - 1)

        # Store wall positions
        self.r_wall = r_wall
        self.t_wall = t_wall
        self.p_wall = p_wall

        # Compute cell centers
        if r_wall[0] == 0.:
            self.r = np.zeros(len(r_wall) - 1)
            self.r[0] = r_wall[1] / 2.
            self.r[1:] = 10. ** ((np.log10(r_wall[1:-1]) + np.log10(r_wall[2:])) / 2.)
        else:
            self.r = 10. ** ((np.log10(r_wall[:-1]) + np.log10(r_wall[1:])) / 2.)

        self.t = (t_wall[:-1] + t_wall[1:]) / 2.
        self.p = (p_wall[:-1] + p_wall[1:]) / 2.

        # Generate 3D versions of r, t, p
        #(each array is 3D and defined in every cell)
        self.gr, self.gt, self.gp = meshgrid_nd(self.r, self.t, self.p)

        # Compute cell centers in cylindrical coordinates
        self.gz = self.gr * np.cos(self.gt)
        self.gw = self.gr * np.sin(self.gt)

        # Generate 3D versions of the inner and outer wall positions respectively
        gr_wall_min, gt_wall_min, gp_wall_min = \
                    meshgrid_nd(r_wall[:-1], t_wall[:-1], p_wall[:-1])
        gr_wall_max, gt_wall_max, gp_wall_max = \
                    meshgrid_nd(r_wall[1:], t_wall[1:], p_wall[1:])

        # USEFUL QUANTITIES

        dr = gr_wall_max - gr_wall_min
        dr2 = gr_wall_max ** 2 - gr_wall_min ** 2
        dr3 = gr_wall_max ** 3 - gr_wall_min ** 3
        dt = gt_wall_max - gt_wall_min
        dcost = np.cos(gt_wall_min) - np.cos(gt_wall_max)
        dp = gp_wall_max - gp_wall_min

        # CELL VOLUMES

        #   dV = dr * (r*dtheta) * (r*sin(theta)*dphi)
        #    V = [r_2^3 - r_1^3] / 3. * [cos(theta_1) - cos(theta_2)] * [phi_2 - phi_1]

        self.volumes = dr3 * dcost * dp / 3.

        # WALL AREAS

        self.areas = np.zeros((6,) + self.shape)

        # R walls:
        #   dA = r^2 * sin(theta) * dtheta * dphi
        #    A = r^2 * [cos(theta_1) - cos(theta_2)] * [phi_2 - phi_1]

        self.areas[0, :, :, :] = gr_wall_min ** 2 * dcost * dp
        self.areas[1, :, :, :] = gr_wall_max ** 2 * dcost * dp

        # Theta walls:
        #   dA = r * sin(theta) * dr * dphi
        #    A = 0.5 * [r_2^2 - r_1^2] * sin(theta) * [phi_2 - phi_1]

        self.areas[2, :, :, :] = 0.5 * dr2 * np.sin(gt_wall_min) * dp
        self.areas[3, :, :, :] = 0.5 * dr2 * np.sin(gt_wall_max) * dp

        # Phi walls:
        #   dA = r * dr * dtheta
        #    A = 0.5 * [r_2^2 - r_1^2] * [theta_2 - theta_1]

        self.areas[4, :, :, :] = 0.5 * dr2 * dt
        self.areas[5, :, :, :] = 0.5 * dr2 * dt

        # CELL WIDTHS

        self.widths = np.zeros((3,) + self.shape)

        # R direction:
        #   dS = dr
        #    S = r_2 - r_1

        self.widths[0, :, :, :] = dr

        # Theta direction:
        #   dS = r * dtheta
        #    S = r * [theta_2 - theta_1]

        self.widths[1, :, :, :] = self.gr * dt

        # Phi direction:
        #   dS = r * sin(theta) * dphi
        #    S = r * sin(theta) * [phi_2 - phi_1]

        self.widths[2, :, :, :] = self.gr * np.sin(self.gt) * dp

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

        if isinstance(array, SphericalPolarGridView):
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
        Read the geometry and physical quantities from a spherical polar grid

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
        Read in geometry information from a spherical polar grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid geometry from.
        '''

        if group.attrs['grid_type'].decode('utf-8') != 'sph_pol':
            raise ValueError("Grid is not spherical polar")

        self.set_walls(group['walls_1']['r'],
                       group['walls_2']['t'],
                       group['walls_3']['p'])

        # Check that advertised hash matches real hash
        if group.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

    def read_quantities(self, group, quantities='all'):
        '''
        Read in physical quantities from a spherical polar grid

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
        Write out the spherical polar grid

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

        g_geometry.attrs['grid_type'] = np.string_('sph_pol'.encode('utf-8'))
        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        dset = g_geometry.create_dataset("walls_1", data=np.array(list(zip(self.r_wall)), dtype=[('r', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = np.string_('cm'.encode('utf-8'))

        dset = g_geometry.create_dataset("walls_2", data=np.array(list(zip(self.t_wall)), dtype=[('t', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = np.string_('rad'.encode('utf-8'))

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
        geo_hash.update(self.r_wall.tostring())
        geo_hash.update(self.t_wall.tostring())
        geo_hash.update(self.p_wall.tostring())
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return SphericalPolarGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, SphericalPolarGridView):
            if self.r_wall is None and self.t_wall is None and self.p_wall is None:
                logger.warn("No geometry in target grid - copying from original grid")
                self.set_walls(value.r_wall, value.t_wall, value.p_wall)
            self.quantities[item] = deepcopy(value.quantities[value.viewed_quantity])
        elif isinstance(value, h5py.ExternalLink):
            self.quantities[item] = value
        elif value == []:
            self.quantities[item] = []
        else:
            raise ValueError('value should be an empty list, and ExternalLink, or a SphericalPolarGridView instance')

    def __contains__(self, item):
        return self.quantities.__contains__(item)

    def reset_quantities(self):
        self.quantities = {}

    def add_derived_quantity(self, name, function):
        if name in self.quantities:
            raise KeyError(name + ' already exists')
        function(self.quantities)


class SphericalPolarGridView(SphericalPolarGrid):

    def __init__(self, grid, quantity):
        self.viewed_quantity = quantity
        SphericalPolarGrid.__init__(self)
        self.set_walls(grid.r_wall, grid.t_wall, grid.p_wall)
        self.quantities = {quantity: grid.quantities[quantity]}

    def append(self, grid):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        grid : 3D Numpy array or SphericalPolarGridView instance
            The grid to copy the quantity from
        '''
        if isinstance(grid, SphericalPolarGridView):
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
            raise ValueError("grid should be a Numpy array or a SphericalPolarGridView instance")

    def add(self, grid):
        '''
        Used to add quantities from another grid

        Parameters
        ----------
        grid : 3D Numpy array or SphericalPolarGridView instance
            The grid to copy the quantity from
        '''
        if type(self.quantities[self.viewed_quantity]) is list:
            raise Exception("need to first specify the item to add to")
        if isinstance(grid, SphericalPolarGridView):
            if type(grid.quantities[grid.viewed_quantity]) is list:
                raise Exception("need to first specify the item to add")
            self._check_array_dimensions(grid.quantities[grid.viewed_quantity])
            self.quantities[self.viewed_quantity] += grid.quantities[grid.viewed_quantity]
        elif isinstance(grid, np.ndarray):
            self._check_array_dimensions(grid)
            self.quantities[self.viewed_quantity] += grid
        else:
            raise ValueError("grid should be a Numpy array or a SphericalPolarGridView instance")

    def __getitem__(self, item):
        if type(item) is int:
            grid = SphericalPolarGridView(self, self.viewed_quantity)
            grid.quantities = {grid.viewed_quantity: grid.quantities[grid.viewed_quantity][item]}
            return grid
        else:
            return SphericalPolarGrid.__getitem__(self, item)

    def __getattr__(self, attribute):
        if attribute == 'array':
            return self.quantities[self.viewed_quantity]
        else:
            return SphericalPolarGrid.__getattr__(self, attribute)
