import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, is_numpy_array, monotonically_increasing


class SphericalPolarGrid(FreezableClass):

    def __init__(self, *args):

        self.nr = None
        self.nt = None
        self.np = None

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

        self.geometry_id = None

        self._freeze()

        if len(args) > 0:
            self.set(*args)

    def set(self, r_wall, t_wall, p_wall):

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
        self.r = 10. ** ((np.log10(r_wall[:-1]) + np.log10(r_wall[1:])) / 2.)
        if r_wall[0] == 0.:
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

        self.areas = np.zeros((6, self.np, self.nt, self.nr))

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

        self.widths = np.zeros((3, self.np, self.nt, self.nr))

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

    def write_physical_array(self, group, array, name, compression=True,
                             physics_dtype=float):
        '''
        Write out a physical quantity defined on a spherical polar grid

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
        Read in a spherical polar grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to read the grid from
        '''

        if group.attrs['grid_type'] != 'sph_pol':
            raise ValueError("Grid is not spherical polar")

        self.geometry_id = group.attrs['geometry']

        self.set(group['Walls 1'], group['Walls 2'], group['Walls 3'])

    def write_geometry(self, group, compression=True, wall_dtype=float):
        '''
        Write out the spherical polar grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to write the grid to
        compression: bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype: type
            The datatype to use to write the wall positions
        '''

        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'sph_pol'

        dset = group.create_dataset("Walls 1", data=np.array(zip(self.r_wall), dtype=[('r', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = group.create_dataset("Walls 2", data=np.array(zip(self.t_wall), dtype=[('t', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'rad'

        dset = group.create_dataset("Walls 3", data=np.array(zip(self.p_wall), dtype=[('p', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'rad'
