import hashlib

import h5py
import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, is_numpy_array, monotonically_increasing, link_or_copy


class CartesianGrid(FreezableClass):

    def __init__(self, *args):

        self.nx = None
        self.ny = None
        self.nz = None

        self.shape = None

        self.x_wall = None
        self.y_wall = None
        self.z_wall = None

        self.x = None
        self.y = None
        self.z = None

        self.gx = None
        self.gy = None
        self.gz = None

        self.volumes = None
        self.areas = None
        self.widths = None

        self.quantities = {}

        self._freeze()

        if len(args) > 0:
            self.set_walls(*args)

    def set_walls(self, x_wall, y_wall, z_wall):

        if type(x_wall) in [list, tuple]:
            x_wall = np.array(x_wall)
        if type(y_wall) in [list, tuple]:
            y_wall = np.array(y_wall)
        if type(z_wall) in [list, tuple]:
            z_wall = np.array(z_wall)

        if not is_numpy_array(x_wall) or x_wall.ndim != 1:
            raise ValueError("x_wall should be a 1-D sequence")
        if not is_numpy_array(y_wall) or y_wall.ndim != 1:
            raise ValueError("y_wall should be a 1-D sequence")
        if not is_numpy_array(z_wall) or z_wall.ndim != 1:
            raise ValueError("z_wall should be a 1-D sequence")

        if not monotonically_increasing(x_wall):
            raise ValueError("x_wall should be monotonically increasing")
        if not monotonically_increasing(y_wall):
            raise ValueError("y_wall should be monotonically increasing")
        if not monotonically_increasing(z_wall):
            raise ValueError("z_wall should be monotonically increasing")

        # Find number of grid cells
        self.nx = len(x_wall) - 1
        self.ny = len(y_wall) - 1
        self.nz = len(z_wall) - 1

        self.shape = (self.nz, self.ny, self.nx)

        # Store wall positions
        self.x_wall = x_wall
        self.y_wall = y_wall
        self.z_wall = z_wall

        # Compute cell centers
        self.x = (x_wall[:-1] + x_wall[1:]) / 2.
        self.y = (y_wall[:-1] + y_wall[1:]) / 2.
        self.z = (z_wall[:-1] + z_wall[1:]) / 2.

        # Generate 3D versions of r, t, p
        #(each array is 3D and defined in every cell)
        self.gx, self.gy, self.gz = meshgrid_nd(self.x, self.y, self.z)

        # Generate 3D versions of the inner and outer wall positions respectively
        gx_wall_min, gy_wall_min, gz_wall_min = \
                    meshgrid_nd(x_wall[:-1], y_wall[:-1], z_wall[:-1])
        gx_wall_max, gy_wall_max, gz_wall_max = \
                    meshgrid_nd(x_wall[1:], y_wall[1:], z_wall[1:])

        # USEFUL QUANTITIES

        dx = gx_wall_max - gx_wall_min
        dy = gy_wall_max - gy_wall_min
        dz = gz_wall_max - gz_wall_min

        # CELL VOLUMES

        self.volumes = dx * dy * dz

        # WALL AREAS

        self.areas = np.zeros((6, self.nz, self.ny, self.nx))

        # X walls:

        self.areas[0, :, :, :] = dy * dz
        self.areas[1, :, :, :] = dy * dz

        # Y walls:

        self.areas[2, :, :, :] = dx * dz
        self.areas[3, :, :, :] = dx * dz

        # Z walls:

        self.areas[4, :, :, :] = dx * dy
        self.areas[5, :, :, :] = dx * dy

        # CELL WIDTHS

        self.widths = np.zeros((3, self.nz, self.ny, self.nx))

        # X direction:

        self.widths[0, :, :, :] = dx

        # Y direction:

        self.widths[1, :, :, :] = dy

        # Z direction:

        self.widths[2, :, :, :] = dz

    def add_quantity(self, name, array, replace=False):
        '''
        Add a physical quantity to the grid

        Parameters
        ----------
        name: str
            The name of the quantity to add
        array: list of np.ndarray or np.ndarray
            The physical quantity to write out. If specified as a list, then
            it is assumed to an array that is defined for different
            populations (e.g. dust or gas).
        replace: bool
            Whether to replace a quantity that is already in the grid
        '''
        self._check_array_dimensions(array)
        if not replace and name in self.quantities:
            raise ValueError("Quantity is already in grid")
        self.quantities[name] = array

    def _check_array_dimensions(self, array):

        if type(array) in [list, tuple]:

            # Check that dimensions are compatible
            for item in array:
                if item.shape != self.shape:
                    raise ValueError("Arrays in list do not have the right "
                                     "dimensions: %s instead of %s"
                                     % (item.shape, self.shape))

            # Convert list of 3D arrays to a single 4D array
            shape = list(self.shape)
            shape.insert(0, len(array))
            array = np.vstack(array).reshape(*shape)

        elif type(array) == np.ndarray:

            if array.shape != self.shape:
                raise ValueError("Array does not have the right "
                                 "dimensions: %s instead of %s"
                                 % (array.shape, self.shape))

        else:

            raise ValueError("Array should be a list or a Numpy array")

    def read(self, group, quantities='all'):
        '''
        Read in a cartesian grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to read the grid from
        quantities: 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Extract HDF5 groups for geometry and physics

        g_geometry = group['Geometry']
        g_physics = group['Physics']

        # Read in geometry

        if g_geometry.attrs['grid_type'] != 'car':
            raise ValueError("Grid is not cartesian")

        self.set(g_geometry['Walls 1'], g_geometry['Walls 2'], g_geometry['Walls 3'])

        # Check that advertised hash matches real hash
        if g_geometry.attrs['geometry'] != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

        # Read in physical quantities

        for quantity in g_physics:
            if quantities == 'all' or quantity in quantities:
                # TODO - if array is 4D, need to convert to list
                self.quantities[quantity] = np.array(g_physics[quantity].array)

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the cartesian grid

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

        if 'Physics' not in group:
            g_physics = group.create_group('Physics')
        else:
            g_physics = group['Physics']

        # Write out geometry

        g_geometry.attrs['grid_type'] = 'car'
        g_geometry.attrs['geometry'] = self.get_geometry_id()

        dset = g_geometry.create_dataset("Walls 1", data=np.array(zip(self.x_wall), dtype=[('x', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = g_geometry.create_dataset("Walls 2", data=np.array(zip(self.y_wall), dtype=[('y', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = g_geometry.create_dataset("Walls 3", data=np.array(zip(self.z_wall), dtype=[('z', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        # Write out physical quantities

        for quantity in self.quantities:
            if quantities == 'all' or quantity in quantities:
                if isinstance(self.quantities[quantity], h5py.ExternalLink):
                    link_or_copy(g_physics, quantity, self.quantities[quantity], copy, absolute_paths=absolute_paths)
                else:
                    self._check_array_dimensions(self.quantities[quantity])
                    dset = g_physics.create_dataset(quantity, data=self.quantities[quantity],
                                                    compression=compression,
                                                    dtype=physics_dtype)
                    dset.attrs['geometry'] = self.get_geometry_id()


    def get_geometry_id(self):
        geo_hash = hashlib.md5()
        geo_hash.update(self.x_wall)
        geo_hash.update(self.y_wall)
        geo_hash.update(self.z_wall)
        return geo_hash.hexdigest()
