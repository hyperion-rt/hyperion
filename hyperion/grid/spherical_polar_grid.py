import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass


def monotonically_increasing(array):
    for i in range(len(array)-1):
        if not array[i+1] > array[i]:
            return False
    return True


class SphericalPolarGrid(FreezableClass):

    def __init__(self, r_wall, t_wall, p_wall):

        # Check that arrays are monotonically increasing
        if not monotonically_increasing(r_wall):
            raise Exception("Array r_wall is not monotonically increasing")

        if not monotonically_increasing(t_wall):
            raise Exception("Array t_wall is not monotonically increasing")

        if not monotonically_increasing(p_wall):
            raise Exception("Array p_wall is not monotonically increasing")

        # Find number of grid cells
        self.nr = len(r_wall) - 1
        self.nt = len(t_wall) - 1
        self.np = len(p_wall) - 1

        self.shape = (self.np, self.nt, self.nr)

        # Store wall positions
        self.r_wall = r_wall
        self.t_wall = t_wall
        self.p_wall = p_wall

        # Compute cell centers
        self.r = 10.**((np.log10(r_wall[:-1]) + np.log10(r_wall[1:])) / 2.)
        if r_wall[0]==0.:
            self.r[0] = r_wall[1] / 2.
        self.t = (t_wall[:-1] + t_wall[1:]) / 2.
        self.p = (p_wall[:-1] + p_wall[1:]) / 2.

        # Generate 3D versions of r, t, p
        #(each array is 3D and defined in every cell)
        self.gr, self.gt, self.gp = meshgrid_nd(self.r, self.t, self.p)

        # Compute cell centers in cylindrical coordinates
        self.gz = self.gr * np.cos(self.gt)
        self.gw = self.gr * np.sin(self.gt)

        # Generate 3D versions of the inner and outer wall positions respectively
        self.gr_wall_min, self.gt_wall_min, self.gp_wall_min = \
                    meshgrid_nd(r_wall[:-1], t_wall[:-1], p_wall[:-1])

        self.gr_wall_max, self.gt_wall_max, self.gp_wall_max = \
                    meshgrid_nd(r_wall[1:], t_wall[1:], p_wall[1:])

        # USEFUL QUANTITIES

        dr = (self.gr_wall_max - self.gr_wall_min)
        dr2 = (self.gr_wall_max**2 - self.gr_wall_min**2)
        dr3 = (self.gr_wall_max**3 - self.gr_wall_min**3)
        dt = self.gt_wall_max - self.gt_wall_min
        dcost = np.cos(self.gt_wall_min) - np.cos(self.gt_wall_max)
        dp = (self.gp_wall_max - self.gp_wall_min)

        # CELL VOLUMES

        #   dV = dr * (r*dtheta) * (r*sin(theta)*dphi)
        #    V = [r_2^3 - r_1^3] / 3. * [cos(theta_1) - cos(theta_2)] * [phi_2 - phi_1]

        self.volumes = dr3 * dcost * dp / 3.

        # WALL AREAS

        self.areas = np.zeros((6, self.np, self.nt, self.nr))

        # R walls:
        #   dA = r^2 * sin(theta) * dtheta * dphi
        #    A = r^2 * [cos(theta_1) - cos(theta_2)] * [phi_2 - phi_1]

        self.areas[0,:,:,:] = self.gr_wall_min**2 * dcost * dp
        self.areas[1,:,:,:] = self.gr_wall_max**2 * dcost * dp

        # Theta walls:
        #   dA = r * sin(theta) * dr * dphi
        #    A = 0.5 * [r_2^2 - r_1^2] * sin(theta) * [phi_2 - phi_1]

        self.areas[2,:,:,:] = 0.5 * dr2 * np.sin(self.gt_wall_min) * dp
        self.areas[3,:,:,:] = 0.5 * dr2 * np.sin(self.gt_wall_max) * dp

        # Phi walls:
        #   dA = r * dr * dtheta
        #    A = 0.5 * [r_2^2 - r_1^2] * [theta_2 - theta_1]

        self.areas[4,:,:,:] = 0.5 * dr2 * dt
        self.areas[5,:,:,:] = 0.5 * dr2 * dt

        # CELL WIDTHS

        self.widths = np.zeros((3, self.np, self.nt, self.nr))

        # R direction:
        #   dS = dr
        #    S = r_2 - r_1

        self.widths[0,:,:,:] = dr

        # Theta direction:
        #   dS = r * dtheta
        #    S = r * [theta_2 - theta_1]

        self.widths[1,:,:,:] = self.gr * dt

        # Phi direction:
        #   dS = r * sin(theta) * dphi
        #    S = r * sin(theta) * [phi_2 - phi_1]

        self.widths[2,:,:,:] = self.gr * np.sin(self.gt) * dp

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
        group.attrs['grid_type'] = 'sph_pol'

        if volumes:
            dset = group.create_dataset("Volumes", data=self.volumes, compression=compression, dtype=geo_dtype)

        if areas:
            dset = group.create_dataset("Areas", data=self.areas, compression=compression, dtype=geo_dtype)

        if widths:
            dset = group.create_dataset("Widths", data=self.widths, compression=compression, dtype=geo_dtype)

        dset = group.create_dataset("Walls 1", data=np.array(zip(self.r_wall), dtype=[('r',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 2", data=np.array(zip(self.t_wall), dtype=[('t',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'rad'

        dset = group.create_dataset("Walls 3", data=np.array(zip(self.p_wall), dtype=[('p',wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'rad'
