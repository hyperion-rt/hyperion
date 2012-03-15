import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, is_numpy_array, monotonically_increasing


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

        self.geometry_id = None

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

    def write_physical_array(self, group, array, name, compression=True,
                             physics_dtype=float):
        '''
        Write out a physical quantity defined on a cartesian grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to write the array to
        array: list of np.ndarray or np.ndarray
            The physical quantity to write out. If specified as a list, then
            it is assumed to an array that is defined for different dust
            types.
        name: str
            The name of the physical array
        compression: bool
            Whether to compress the array in the HDF5 file
        physics_dtype: type
            The datatype to use to write the array
        '''

        if type(array) in [list, tuple]:

            # Check that dimensions are compatible
            for item in array:
                if item.shape != self.shape:
                    raise ValueError("Grids in list do not have the right "
                                     "dimensions: %s instead of %s"
                                     % (item.shape, self.shape))

            # Convert list of 3D arrays to a single 4D array
            shape = list(self.shape)
            shape.insert(0, len(array))
            array = np.vstack(array).reshape(*shape)

        elif type(array) == np.ndarray:

            if array.shape != self.shape:
                raise ValueError("Grid does not have the right "
                                 "dimensions: %s instead of %s"
                                 % (array.shape, self.shape))

        else:

            raise ValueError("array should be a list or a Numpy array")

        dset = group.create_dataset(name, data=array,
                                    compression=compression,
                                    dtype=physics_dtype)

        dset.attrs['geometry'] = self.geometry_id

    def read_geometry(self, group):
        '''
        Read in a cartesian grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to read the grid from
        '''

        if group.attrs['grid_type'] != 'car':
            raise ValueError("Grid is not cartesian")

        self.geometry_id = group.attrs['geometry']

        self.set(group['Walls 1'], group['Walls 2'], group['Walls 3'])

    def write_geometry(self, group, compression=True, wall_dtype=float):
        '''
        Write out the cartesian grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to write the grid to
        compression: bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype: type
            The datatype to use to write the wall positions
        '''

        group.attrs['grid_type'] = 'car'
        group.attrs['geometry'] = self.geometry_id

        dset = group.create_dataset("Walls 1", data=np.array(zip(self.x_wall), dtype=[('x', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 2", data=np.array(zip(self.y_wall), dtype=[('y', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 3", data=np.array(zip(self.z_wall), dtype=[('z', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'
