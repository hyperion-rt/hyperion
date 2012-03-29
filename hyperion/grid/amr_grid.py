import struct
import hashlib

import h5py
import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, link_or_copy


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
            return FreezableClass.__getattr__(self, attribute)


class Level(FreezableClass):

    def __init__(self):

        # The list of grids in the level
        self.grids = []

        self._freeze()


class AMRGrid(FreezableClass):

    def __init__(self):

        self.levels = []

        self.geometry_id = None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (1, 1, self.ncells)
        else:
            return FreezableClass(self, attribute)

    def _check_array_dimensions(self, amr_grid=None):
        '''
        Check that a grid's array dimensions agree with this grid's metadata

        Parameters
        ----------
        amr_grid: AMR grid, optional
            The AMR grid for which to test the array dimensions. If this is not
            specified, this method performs a self-consistency check of array
            dimensions and meta-data.
        '''

        # If no grid is specified, do a self-consistency checks
        if amr_grid is None:
            amr_grid = self

        # Loop over levels
        for ilevel, level_ref in enumerate(self.levels):

            # Read in level
            level = amr_grid.levels[ilevel]

            # Loop over grids
            for igrid, grid_ref in enumerate(level.grids):

                # Read in grid
                grid = level.grids[igrid]

                # Loop over quantities
                for quantity in grid.quantities:

                    # Extract array
                    array = grid.quantities[quantity]

                    if type(array) in [list, tuple]:

                        # Check that dimensions are compatible
                        for item in array:
                            if item.shape != grid_ref.shape:
                                raise ValueError("Arrays in list do not have the right "
                                                 "dimensions: %s instead of %s"
                                                 % (item.shape, grid_ref.shape))

                    elif type(array) == np.ndarray:

                        if array.shape != grid_ref.shape:
                            raise ValueError("Array does not have the right "
                                             "dimensions: %s instead of %s"
                                             % (array.shape, grid_ref.shape))

                    else:

                        raise ValueError("Array should be a list or a Numpy array")

    def read(self, group, quantities='all'):
        '''
        Read in an AMR grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to read the grid from
        quantities: 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Initialize levels
        self.levels = []

        # Read geometry group
        g_geometry = group['Geometry']

        # Read quantities group
        g_quantities = group['Quantities']

        # Read geometry ID
        self.geometry_id = g_geometry.attrs['geometry']

        # Check that grid is indeed AMR
        if g_geometry.attrs['grid_type'] != 'amr':
            raise Exception("Grid is not AMR")

        # Initialize levels list
        self.levels = []

        # Loop over levels
        for ilevel in range(g_geometry.attrs['nlevels']):

            # Read in level
            level_path = 'level_%05i' % (ilevel + 1)
            g_level = g_geometry[level_path]

            # Initialize level
            level = Level()

            # Loop over grids
            for igrid in range(g_level.attrs['ngrids']):

                # Read in grid
                grid_path = 'grid_%05i' % (igrid + 1)
                g_grid = g_level[grid_path]

                # Initialize grid
                grid = Grid()

                # Retrieve real-world grid boundaries
                grid.xmin = g_grid.attrs['xmin']
                grid.xmax = g_grid.attrs['xmax']
                grid.ymin = g_grid.attrs['ymin']
                grid.ymax = g_grid.attrs['ymax']
                grid.zmin = g_grid.attrs['zmin']
                grid.zmax = g_grid.attrs['zmax']

                # Retrieve grid dimensions
                grid.nx = g_grid.attrs['n1']
                grid.ny = g_grid.attrs['n2']
                grid.nz = g_grid.attrs['n3']

                # Read in desired quantities
                g_grid_quantities = g_quantities[level_path][grid_path]
                for quantity in g_grid_quantities:
                    if quantities == 'all' or quantity in quantities:
                        # TODO - if array is 4D, need to convert to list
                        grid.quantities[quantity] = np.array(g_grid_quantities[quantity])

                # Append grid to current level
                level.grids.append(grid)

            # Append level to overall grid
            self.levels.append(level)

        # Check that advertised hash matches real hash
        if g_geometry.attrs['geometry'] != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the AMR grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to write the grid to
        quantities: 'all' or list
            Which physical quantities to write out. Use 'all' to write out all
            quantities or a list of strings to write only specific quantities.
        compression: bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype: type
            The datatype to use to write the wall positions
        physics_dtype: type
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

        g_geometry.attrs['grid_type'] = 'amr'
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

        g_geometry.attrs['geometry'] = self.get_geometry_id()

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
