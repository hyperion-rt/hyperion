import numpy as np

from hyperion.util.functions import FreezableClass


class OcTreeGrid(FreezableClass):

    def __init__(self, refined, x, y, z, dx, dy, dz):

        self.refined = refined
        # Find number of grid cells
        self.ncells = len(self.refined)

        self.x = x
        self.y = y
        self.z = z
        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.geometry_id = None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (1, 1, self.ncells)

    def write_physical_array(self, group, array, name, dust=False, compression=True, physics_dtype=float):

        if dust:
            shape = list(self.shape)
            shape.insert(0, len(array))
            array = np.vstack(array).reshape(*shape)
        dset = group.create_dataset(name, data=array, compression=compression, dtype=physics_dtype)
        dset.attrs['geometry'] = self.geometry_id

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):

        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'oct'

        # if volumes:
        #     dset = group.create_dataset("Volumes", data=self.volumes, compression=compression, dtype=geo_dtype)
        #
        # if areas:
        #     dset = group.create_dataset("Areas", data=self.areas, compression=compression, dtype=geo_dtype)
        #
        # if widths:
        #     dset = group.create_dataset("Widths", data=self.widths, compression=compression, dtype=geo_dtype)

        group.attrs['x'] = self.x
        group.attrs['y'] = self.y
        group.attrs['z'] = self.z
        group.attrs['dx'] = self.dx
        group.attrs['dy'] = self.dy
        group.attrs['dz'] = self.dz

        dset = group.create_dataset("Cells", data=np.array(zip(self.refined), dtype=[('refined', np.int32)]), compression=compression)
