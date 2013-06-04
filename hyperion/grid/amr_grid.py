from __future__ import print_function, division

import os
import struct
import hashlib
from copy import deepcopy

import h5py
import numpy as np

from ..util.meshgrid import meshgrid_nd
from ..util.functions import FreezableClass, link_or_copy
from astropy import log as logger
from .grid_helpers import single_grid_dims


def zero_density(grid, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=np.inf, zmax=np.inf):
    for ilevel, level in enumerate(grid.levels):
        for igrid, grid in enumerate(level.grids):
            wx = np.linspace(grid.xmin, grid.xmax, grid.nx + 1)
            wy = np.linspace(grid.ymin, grid.ymax, grid.ny + 1)
            wz = np.linspace(grid.zmin, grid.zmax, grid.nz + 1)
            x = 0.5 * (wx[:-1] + wx[1:])
            y = 0.5 * (wy[:-1] + wy[1:])
            z = 0.5 * (wz[:-1] + wz[1:])
            gx, gy, gz = meshgrid_nd(x, y, z)
            reset = (gx < xmin) | (gx > xmax) | (gy < ymin) | (gy > ymax) | (gz < zmin) | (gz > zmax)
            grid.data[reset] = 0.
    return grid


class Grid(FreezableClass):

    def __init__(self):

        # The 3D data arrays (can have various components)
        self.quantities = {}

        # The boundaries of the 3D grid in real space
        self.xmin, self.xmax = None, None
        self.ymin, self.ymax = None, None
        self.zmin, self.zmax = None, None

        # The dimensions of the array
        self.nx, self.ny, self.nz = None, None, None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (self.nz, self.ny, self.nx)
        else:
            return FreezableClass.__getattribute__(self, attribute)


class Level(FreezableClass):

    def __init__(self):

        # The list of grids in the level
        self.grids = []

        self._freeze()

    def add_grid(self):
        grid = Grid()
        self.grids.append(grid)
        return grid


class AMRGrid(FreezableClass):
    '''
    An AMR grid.

    Levels are stored in the ``levels`` attribute, which is a list of
    :class:`hyperion.grid.amr_grid.Level` objects, which in turn
    contain a ``grids`` attribute which is a list of
    :class:`~hyperion.grid.amr_grid.Grid` objects.

    Levels can be added with::

        level = amr.add_level()

    And grids can be added to a level with::

        grid = level.add_grid()

    Grid objects have the following attributes which should be set:

        * ``xmin`` - lower x position of the grid
        * ``xmax`` - upper x position of the grid
        * ``ymin`` - lower y position of the grid
        * ``ymax`` - upper y position of the grid
        * ``zmin`` - lower z position of the grid
        * ``zmax`` - upper z position of the grid
        * ``nx`` - number of cells in x direction
        * ``ny`` - number of cells in y direction
        * ``nz`` - number of cells in z direction
        * ``quantities`` - a dictionary containing physical quantities (see below)

    :class:`~hyperion.grid.AMRGrid` objects may contain multiple
    quantities (e.g. density, specific energy). To access these, you can
    specify the name of the quantity as an item::

         >>> grid['density']

    which is no longer an :class:`~hyperion.grid.AMRGrid` object, but
    a :class:`~hyperion.grid.AMRGridView` object. When setting
    this for the first time, this can be set either to another
    :class:`~hyperion.grid.AMRGridView` object, an external h5py
    link, or an empty list. For example, the following should work:

        >>> grid['density_new'] = grid['density']

    :class:`~hyperion.grid.AMRGridView` objects allow the
    specific dust population to be selected as an index:

        >>> grid['density'][0]

    Which is also an :class:`~hyperion.grid.AMRGridView` object.
    '''

    def __init__(self, amr_grid=None):

        # Initalize AMR levels
        self.levels = []

        self._freeze()

        # Copy geometry if provided
        if amr_grid is not None:
            for level in amr_grid.levels:
                level_ref = self.add_level()
                for grid in level.grids:
                    grid_ref = level_ref.add_grid()
                    grid_ref.nx = grid.nx
                    grid_ref.ny = grid.ny
                    grid_ref.nz = grid.nz
                    grid_ref.xmin, grid_ref.xmax = grid.xmin, grid.xmax
                    grid_ref.ymin, grid_ref.ymax = grid.ymin, grid.ymax
                    grid_ref.zmin, grid_ref.zmax = grid.zmin, grid.zmax
                    grid_ref.quantities = {}

    def remove_level(self, level_id):
        self.levels.pop(level_id)

    def add_level(self):
        level = Level()
        self.levels.append(level)
        return level

    def __getattr__(self, attribute):
        if attribute == 'n_dust':
            n_dust = None
            for level in self.levels:
                for grid in level.grids:
                    for quantity in grid.quantities:
                        n_dust_q, shape_q = single_grid_dims(grid.quantities[quantity])
                        if n_dust is None:
                            n_dust = n_dust_q
                        elif n_dust_q is not None:
                            if n_dust != n_dust_q:
                                raise ValueError("Not all dust lists in the grid have the same size")
            return n_dust
        else:
            return FreezableClass.__getattribute__(self, attribute)

    def _check_array_dimensions(self, amr_grid=None):
        '''
        Check that a grid's array dimensions agree with this grid's metadata

        Parameters
        ----------
        amr_grid : AMR grid, optional
            The AMR grid for which to test the array dimensions. If this is not
            specified, this method performs a self-consistency check of array
            dimensions and meta-data.
        '''

        # If no grid is specified, do a self-consistency check
        if amr_grid is None:
            amr_grid = self

        n_pop_ref = None

        # Loop over levels
        for ilevel, level_ref in enumerate(self.levels):

            # Read in level
            level = amr_grid.levels[ilevel]

            # Loop over grids
            for igrid, grid_ref in enumerate(level_ref.grids):

                # Read in grid
                grid = level.grids[igrid]

                # Loop over quantities
                for quantity in grid.quantities:

                    n_pop, shape = single_grid_dims(grid.quantities[quantity])

                    if shape != grid_ref.shape:
                        raise ValueError("Quantity arrays do not have the right "
                                         "dimensions: %s instead of %s"
                                         % (shape, grid_ref.shape))

                    if n_pop is not None:
                        if n_pop_ref is None:
                            n_pop_ref = n_pop
                        elif n_pop != n_pop_ref:
                            raise ValueError("Not all dust lists in the grid have the same size")

    def read(self, group, quantities='all'):
        '''
        Read the geometry and physical quantities from an AMR grid

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
        Read in geometry information from an AMR grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the geometry from
        '''

        # Check that grid is indeed AMR
        if group.attrs['grid_type'].decode('utf-8') != 'amr':
            raise Exception("Grid is not an AMR grid")

        # Initialize levels list
        self.levels = []

        # Loop over levels
        for ilevel in range(group.attrs['nlevels']):

            # Read in level
            level_path = 'level_%05i' % (ilevel + 1)
            g_level = group[level_path]

            # Initialize level
            level = self.add_level()

            # Loop over grids
            for igrid in range(g_level.attrs['ngrids']):

                # Read in grid
                grid_path = 'grid_%05i' % (igrid + 1)
                g_grid = g_level[grid_path]

                # Initialize grid
                grid = level.add_grid()

                # Retrieve real-world grid boundaries
                grid.xmin = g_grid.attrs['xmin']
                grid.xmax = g_grid.attrs['xmax']
                grid.ymin = g_grid.attrs['ymin']
                grid.ymax = g_grid.attrs['ymax']
                grid.zmin = g_grid.attrs['zmin']
                grid.zmax = g_grid.attrs['zmax']

                # Retrieve grid dimensions
                grid.nx = int(g_grid.attrs['n1'])
                grid.ny = int(g_grid.attrs['n2'])
                grid.nz = int(g_grid.attrs['n3'])

        # Check that advertised hash matches real hash
        if group.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

    def read_quantities(self, group, quantities='all'):
        '''
        Read in physical quantities from an AMR grid

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to read the grid quantities from
        quantities : 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Loop over levels
        for ilevel, level in enumerate(self.levels):

            # Read in level
            level_path = 'level_%05i' % (ilevel + 1)

            # Loop over grids
            for igrid, grid in enumerate(level.grids):

                # Read in grid
                grid_path = 'grid_%05i' % (igrid + 1)

                # Read in desired quantities
                g_grid_quantities = group[level_path][grid_path]
                for quantity in g_grid_quantities:
                    if quantities == 'all' or quantity in quantities:
                        array = np.array(g_grid_quantities[quantity])
                        if array.ndim == 4:  # if array is 4D, it is a list of 3D arrays
                            grid.quantities[quantity] = [array[i] for i in range(array.shape[0])]
                        else:
                            grid.quantities[quantity] = array

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the AMR grid

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

        g_geometry.attrs['grid_type'] = np.string_('amr'.encode('utf-8'))
        g_geometry.attrs['nlevels'] = len(self.levels)

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

        # Write out physical quantities

        # Loop over levels
        for ilevel, level in enumerate(self.levels):

            # Read in level
            level_path = 'level_%05i' % (ilevel + 1)
            g_level = g_geometry.create_group(level_path)
            q_level = g_quantities.create_group(level_path)
            g_level.attrs['ngrids'] = len(level.grids)

            # Loop over grids
            for igrid, grid in enumerate(level.grids):

                # Read in grid
                grid_path = 'grid_%05i' % (igrid + 1)
                g_grid = g_level.create_group(grid_path)
                q_grid = q_level.create_group(grid_path)

                # Write real-world grid boundaries
                g_grid.attrs['xmin'] = grid.xmin
                g_grid.attrs['xmax'] = grid.xmax
                g_grid.attrs['ymin'] = grid.ymin
                g_grid.attrs['ymax'] = grid.ymax
                g_grid.attrs['zmin'] = grid.zmin
                g_grid.attrs['zmax'] = grid.zmax

                # Wrote grid dimensions
                g_grid.attrs['n1'] = grid.nx
                g_grid.attrs['n2'] = grid.ny
                g_grid.attrs['n3'] = grid.nz

                # Write out physical quantities
                for quantity in grid.quantities:
                    if quantities == 'all' or quantity in quantities:
                        if isinstance(grid.quantities[quantity], h5py.ExternalLink):
                            link_or_copy(q_grid, quantity, grid.quantities[quantity], copy, absolute_paths=absolute_paths)
                        else:
                            q_grid.create_dataset(quantity, data=grid.quantities[quantity],
                                                  compression=compression,
                                                  dtype=physics_dtype)

        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

    def write_single_array(self, group, name, amr_grid, copy=True, absolute_paths=False, compression=True, physics_dtype=float):
        '''
        Write out a single quantity, checking for consistency with geometry

        Parameters
        ----------
        group : h5py.Group
            The HDF5 group to write the grid to
        name : str
            The name of the array in the group
        amr_grid : AMRGridView
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

        if not isinstance(amr_grid, AMRGridView):
            raise ValueError("amr_grid should be an AMRGridView instance")

        # Loop over levels
        for ilevel, level in enumerate(self.levels):

            # Read in level
            level_path = 'level_%05i' % (ilevel + 1)
            if level_path in group:
                q_level = group[level_path]
            else:
                q_level = group.create_group(level_path)

            # Loop over grids
            for igrid, grid in enumerate(level.grids):

                # Read in grid
                grid_path = 'grid_%05i' % (igrid + 1)
                if grid_path in q_level:
                    q_grid = q_level[grid_path]
                else:
                    q_grid = q_level.create_group(grid_path)

                # Write out physical quantities
                array = amr_grid.levels[ilevel].grids[igrid].quantities[amr_grid.viewed_quantity]
                if isinstance(array, h5py.ExternalLink):
                    link_or_copy(q_grid, name, array, copy, absolute_paths=absolute_paths)
                else:
                    q_grid.create_dataset(name, data=array,
                                          compression=compression,
                                          dtype=physics_dtype)

    def get_geometry_id(self):
        geo_hash = hashlib.md5()
        for level in self.levels:
            for grid in level.grids:
                geo_hash.update(struct.pack('>d', grid.xmin))
                geo_hash.update(struct.pack('>d', grid.xmax))
                geo_hash.update(struct.pack('>d', grid.ymin))
                geo_hash.update(struct.pack('>d', grid.ymax))
                geo_hash.update(struct.pack('>d', grid.zmin))
                geo_hash.update(struct.pack('>d', grid.zmax))
                geo_hash.update(struct.pack('>q', grid.nx))
                geo_hash.update(struct.pack('>q', grid.ny))
                geo_hash.update(struct.pack('>q', grid.nz))
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return AMRGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, AMRGridView):
            if self.levels == [] and value.levels != []:
                logger.warn("No geometry in target grid - copying from original grid")
                for level in value.levels:
                    level_ref = self.add_level()
                    for grid in level.grids:
                        grid_ref = level_ref.add_grid()
                        grid_ref.nx = grid.nx
                        grid_ref.ny = grid.ny
                        grid_ref.nz = grid.nz
                        grid_ref.xmin, grid_ref.xmax = grid.xmin, grid.xmax
                        grid_ref.ymin, grid_ref.ymax = grid.ymin, grid.ymax
                        grid_ref.zmin, grid_ref.zmax = grid.zmin, grid.zmax
                        grid_ref.quantities = {}
            for ilevel, level_ref in enumerate(self.levels):
                level = value.levels[ilevel]
                for igrid, grid_ref in enumerate(level_ref.grids):
                    grid = level.grids[igrid]
                    grid_ref.quantities[item] = deepcopy(grid.quantities[value.viewed_quantity])
        elif isinstance(value, h5py.ExternalLink):
            filename = value.filename
            base_path = os.path.dirname(value.path)
            array_name = os.path.basename(value.path)
            for ilevel, level_ref in enumerate(self.levels):
                level_path = 'level_%05i' % (ilevel + 1)
                for igrid, grid_ref in enumerate(level_ref.grids):
                    grid_path = 'grid_%05i' % (ilevel + 1)
                    grid_ref.quantities[item] = h5py.ExternalLink(filename, os.path.join(base_path, level_path, grid_path, array_name))
        elif value == []:
            for level in self.levels:
                for grid in level.grids:
                    grid.quantities[item] = []
        else:
            raise ValueError('value should be an empty list or an AMRGridView instance')

    def __contains__(self, item):
        if len(self.levels) > 0:
            if len(self.levels[0].grids) > 0:
                return item in self.levels[0].grids[0].quantities
            else:
                return False
        else:
            return False

    def reset_quantities(self):
        self.quantities = {}
        for level in self.levels:
            for grid in level.grids:
                grid.quantities = {}

    def add_derived_quantity(self, name, function):
        for level in self.levels:
            for grid in level.grids:
                if name in grid.quantities:
                    raise KeyError(name + ' already exists')
                function(grid.quantities)

    def to_yt(self, dust_id=0):
        '''
        Convert AMR grid to a yt object (requires yt)

        Parameters
        ----------
        dust_id : int, optional
            The ID of the dust population to extract. If not set, this
            defaults to 0 (the first dust population).
        '''
        from yt_wrappers import amr_grid_to_yt_stream
        return amr_grid_to_yt_stream(self.levels, dust_id)


class AMRGridView(AMRGrid):

    def __init__(self, amr_grid, quantity):

        self.viewed_quantity = quantity
        AMRGrid.__init__(self)

        for level_ref in amr_grid.levels:
            level = self.add_level()
            for grid_ref in level_ref.grids:
                grid = level.add_grid()
                grid.nx = grid_ref.nx
                grid.ny = grid_ref.ny
                grid.nz = grid_ref.nz
                grid.xmin, grid.xmax = grid_ref.xmin, grid_ref.xmax
                grid.ymin, grid.ymax = grid_ref.ymin, grid_ref.ymax
                grid.zmin, grid.zmax = grid_ref.zmin, grid_ref.zmax
                grid.quantities = {}
                grid.quantities[quantity] = grid_ref.quantities[quantity]

    def append(self, amr_grid_view):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        amr_grid : AMRGridView instance
            The grid to copy the quantity from
        '''
        if not isinstance(amr_grid_view, AMRGridView):
            raise ValueError("amr_grid_view should be an AMRGridView instance")
        self._check_array_dimensions(amr_grid_view[amr_grid_view.viewed_quantity])
        for ilevel, level_ref in enumerate(self.levels):
            level = amr_grid_view.levels[ilevel]
            for igrid, grid_ref in enumerate(level_ref.grids):
                grid = level.grids[igrid]
                if grid_ref.quantities[self.viewed_quantity] is grid.quantities[amr_grid_view.viewed_quantity]:
                    raise Exception("Calling append recursively")
                if type(grid.quantities[amr_grid_view.viewed_quantity]) is list:
                    raise Exception("Can only append a single grid")
                grid_ref.quantities[self.viewed_quantity].append(deepcopy(grid.quantities[amr_grid_view.viewed_quantity]))

    def add(self, amr_grid_view):
        '''
        Used to add quantities from another grid

        Parameters
        ----------
        amr_grid : AMRGridView instance
            The grid to copy the quantity from
        '''
        if not isinstance(amr_grid_view, AMRGridView):
            raise ValueError("amr_grid_view should be an AMRGridView instance")
        self._check_array_dimensions(amr_grid_view[amr_grid_view.viewed_quantity])
        for ilevel, level_ref in enumerate(self.levels):
            level = amr_grid_view.levels[ilevel]
            for igrid, grid_ref in enumerate(level_ref.grids):
                grid = level.grids[igrid]
                if type(grid_ref.quantities[self.viewed_quantity]) is list:
                    raise Exception("need to first specify the item to add to")
                if type(grid.quantities[amr_grid_view.viewed_quantity]) is list:
                    raise Exception("need to first specify the item to add")
                grid_ref.quantities[self.viewed_quantity] += grid.quantities[amr_grid_view.viewed_quantity]

    def __getitem__(self, item):
        if type(item) is int:
            amr_grid = AMRGridView(self, self.viewed_quantity)
            for level in amr_grid.levels:
                for grid in level.grids:
                    grid.quantities = {amr_grid.viewed_quantity: grid.quantities[amr_grid.viewed_quantity][item]}
            return amr_grid
        else:
            return AMRGrid.__getitem__(self, item)
