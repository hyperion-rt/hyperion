import hashlib

import h5py
import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass, is_numpy_array, monotonically_increasing, link_or_copy


class CylindricalPolarGrid(FreezableClass):

    def __init__(self, *args):

        self.nw = None
        self.nz = None
        self.np = None

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

        self.areas = np.zeros((6, self.np, self.nz, self.nw))

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
        Read in a cylindrical polar grid

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

        if g_geometry.attrs['grid_type'] != 'cyl_pol':
            raise ValueError("Grid is not cylindrical polar")

        self.set(g_geometry['Walls 1'], g_geometry['Walls 2'], g_geometry['Walls 3'])

        # Check that advertised hash matches real hash
        if g_geometry.attrs['geometry'] != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

        # Read in physical quantities

        for quantity in g_physics:
            if quantities == 'all' or quantity in quantities:
                self.quantities[quantity] = np.array(g_physics[quantity].array)

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the cylindrical polar grid

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

        g_geometry.attrs['grid_type'] = 'cyl_pol'
        g_geometry.attrs['geometry'] = self.get_geometry_id()

        dset = g_geometry.create_dataset("Walls 1", data=np.array(zip(self.w_wall), dtype=[('w', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = g_geometry.create_dataset("Walls 2", data=np.array(zip(self.z_wall), dtype=[('z', wall_dtype)]), compression=compression)
        dset.attrs['Unit'] = 'cm'

        dset = g_geometry.create_dataset("Walls 3", data=np.array(zip(self.p_wall), dtype=[('p', wall_dtype)]), compression=compression)
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
        geo_hash.update(self.w_wall)
        geo_hash.update(self.z_wall)
        geo_hash.update(self.p_wall)
        return geo_hash.hexdigest()
