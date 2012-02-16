import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, is_numpy_array, monotonically_increasing


class CylindricalPolarGrid(FreezableClass):

    def __init__(self, w_wall, z_wall, p_wall):

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

        # Find number of grid cells
        self.nw = len(w_wall) - 1
        self.nz = len(z_wall) - 1
        self.np = len(p_wall) - 1

        self.shape = (self.np, self.nz, self.nw)

        # Store wall positions
        self.w_wall = w_wall
        self.z_wall = z_wall
        self.p_wall = p_wall

        # Compute cell centers
        self.w = 10. ** ((np.log10(w_wall[:-1]) + np.log10(w_wall[1:])) / 2.)
        if w_wall[0] == 0.:
            self.w[0] = w_wall[1] / 2.
        self.z = (z_wall[:-1] + z_wall[1:]) / 2.
        self.p = (p_wall[:-1] + p_wall[1:]) / 2.

        # Generate 3D versions of r, t, p
        #(each array is 3D and defined in every cell)
        self.gw, self.gz, self.gp = meshgrid_nd(self.w, self.z, self.p)

        # Generate 3D versions of the inner and outer wall positions respectively
        self.gw_wall_min, self.gz_wall_min, self.gp_wall_min = \
                    meshgrid_nd(w_wall[:-1], z_wall[:-1], p_wall[:-1])

        self.gw_wall_max, self.gz_wall_max, self.gp_wall_max = \
                    meshgrid_nd(w_wall[1:], z_wall[1:], p_wall[1:])

        # USEFUL QUANTITIES

        dr = self.gw_wall_max - self.gw_wall_min
        dr2 = (self.gw_wall_max ** 2 - self.gw_wall_min ** 2)
        dz = self.gz_wall_max - self.gz_wall_min
        dp = self.gp_wall_max - self.gp_wall_min

        # CELL VOLUMES

        #   dV = dr * dz * (r*dphi)
        #    V = [r_2^2 - r_1^2] / 2. * [z_2 - z_1] * [phi_2 - phi_1]

        self.volumes = dr2 * dz * dp / 2.

        # WALL AREAS

        self.areas = np.zeros((6, self.np, self.nz, self.nw))

        # R walls:
        #   dA = r * dz * dphi
        #    A = r * [z 2 - z_1] * [phi_2 - phi_1]

        self.areas[0, :, :, :] = self.gw_wall_min * dz * dp
        self.areas[1, :, :, :] = self.gw_wall_max * dz * dp

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

        self.widths = np.zeros((3, self.np, self.nz, self.nw))

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

        self.geometry_id = None

        self._freeze()

    def write_physical_array(self, group, array, name, dust=False, compression=True, physics_dtype=float):

        if dust:
            shape = list(self.shape)
            shape.insert(0, len(array))
            array = np.vstack(array).reshape(*shape)
        dset = group.create_dataset(name, data=array, compression=compression, dtype=physics_dtype)
        dset.attrs['geometry'] = self.geometry_id

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):

        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'cyl_pol'

        if volumes:
            dset = group.create_dataset("Volumes", data=self.volumes, compression=compression, dtype=geo_dtype)

        if areas:
            dset = group.create_dataset("Areas", data=self.areas, compression=compression, dtype=geo_dtype)

        if widths:
            dset = group.create_dataset("Widths", data=self.widths, compression=compression, dtype=geo_dtype)

        dset = group.create_dataset("Walls 1", data=np.array(zip(self.w_wall), dtype=[('w', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 2", data=np.array(zip(self.z_wall), dtype=[('z', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 3", data=np.array(zip(self.p_wall), dtype=[('p', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'rad'
