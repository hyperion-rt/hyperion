import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass


class CartesianGrid(FreezableClass):

    def __init__(self, x_wall, y_wall, z_wall):

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
        self.gx_wall_min, self.gy_wall_min, self.gz_wall_min = \
                    meshgrid_nd(x_wall[:-1], y_wall[:-1], z_wall[:-1])

        self.gx_wall_max, self.gy_wall_max, self.gz_wall_max = \
                    meshgrid_nd(x_wall[1:], y_wall[1:], z_wall[1:])

        # USEFUL QUANTITIES

        dx = (self.gx_wall_max - self.gx_wall_min)
        dy = (self.gy_wall_max - self.gy_wall_min)
        dz = (self.gz_wall_max - self.gz_wall_min)

        # CELL VOLUMES

        self.volumes = dx * dy * dz

        # WALL AREAS

        self.areas = np.zeros((6, self.nz, self.ny, self.nx))

        # X walls:

        self.areas[0,:,:,:] = dy * dz
        self.areas[1,:,:,:] = dy * dz

        # Y walls:

        self.areas[2,:,:,:] = dx * dz
        self.areas[3,:,:,:] = dx * dz

        # Z walls:

        self.areas[4,:,:,:] = dx * dy
        self.areas[5,:,:,:] = dx * dy

        # CELL WIDTHS

        self.widths = np.zeros((3, self.nz, self.ny, self.nx))

        # X direction:

        self.widths[0,:,:,:] = dx

        # Y direction:

        self.widths[1,:,:,:] = dy

        # Z direction:

        self.widths[2,:,:,:] = dz

        self.geometry_id = None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (self.nx, self.ny, self.nz)

    def write_physical_array(self, group, array, name, dust=False, compression=True, physics_dtype=float):

        if dust:
            shape = list(self.shape)
            shape.insert(0, len(array))
            array = np.vstack(array).reshape(*shape)
        dset = group.create_dataset(name, data=array, compression=compression, dtype=physics_dtype)
        dset.attrs['geometry'] = self.geometry_id

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):

        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'car'

        if volumes:
            dset = group.create_dataset("Volumes", data=self.volumes, compression=compression, dtype=geo_dtype)

        if areas:
            dset = group.create_dataset("Areas", data=self.areas, compression=compression, dtype=geo_dtype)

        if widths:
            dset = group.create_dataset("Widths", data=self.widths, compression=compression, dtype=geo_dtype)

        dset = group.create_dataset("Walls 1", data=np.array(zip(self.x_wall), dtype=[('x',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 2", data=np.array(zip(self.y_wall), dtype=[('y',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 3", data=np.array(zip(self.z_wall), dtype=[('z',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'
