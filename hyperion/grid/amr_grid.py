import numpy as np

from hyperion.util.functions import FreezableClass


class AMRGrid(FreezableClass):

    def __init__(self, description):

        self.levels = description

        self.geometry_id = None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (1, 1, self.ncells)

    def write_physical_array(self, group, array, name, dust=False, compression=True, physics_dtype=float):

        for ilevel,level in enumerate(self.levels):
            level_name = "Level %i" % (ilevel+1)
            if level_name in group:
                g_level = group[level_name]
            else:
                g_level = group.create_group("Level %i" % (ilevel+1))
            for ifab, fab in enumerate(level):
                fab_name = "Fab %i" % (ifab+1)
                if fab_name in g_level:
                    g_fab = g_level[fab_name]
                else:
                    g_fab = g_level.create_group("Fab %i" % (ifab+1))
                if dust:
                    shape = list(array[0][ilevel][ifab].shape)
                    shape.insert(0, len(array))
                    fab_array = np.vstack([a[ilevel][ifab] for a in array]).reshape(*shape)
                else:
                    fab_array = array[ilevel][ifab]
                g_fab.create_dataset(name, data=fab_array, compression=compression, dtype=physics_dtype)

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):

        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'amr'

        # if volumes:
        #     dset = group.create_dataset("Volumes", data=self.volumes, compression=compression, dtype=geo_dtype)
        #
        # if areas:
        #     dset = group.create_dataset("Areas", data=self.areas, compression=compression, dtype=geo_dtype)
        #
        # if widths:
        #     dset = group.create_dataset("Widths", data=self.widths, compression=compression, dtype=geo_dtype)

        group.attrs['nlevels'] = len(self.levels)
        for ilevel, level in enumerate(self.levels):

            wx = []
            wy = []
            wz = []

            g_level = group.create_group("Level %i" % (ilevel+1))
            g_level.attrs['nfabs'] = len(level)
            for ifab, fab in enumerate(level):
                g_fab = g_level.create_group("Fab %i" % (ifab+1))
                g_fab.attrs['xmin'] = fab[0]
                g_fab.attrs['xmax'] = fab[1]
                g_fab.attrs['ymin'] = fab[2]
                g_fab.attrs['ymax'] = fab[3]
                g_fab.attrs['zmin'] = fab[4]
                g_fab.attrs['zmax'] = fab[5]
                g_fab.attrs['n1'] = fab[6]
                g_fab.attrs['n2'] = fab[7]
                g_fab.attrs['n3'] = fab[8]
                wx.append(np.linspace(fab[0], fab[1], fab[6]+1))
                wy.append(np.linspace(fab[2], fab[3], fab[7]+1))
                wz.append(np.linspace(fab[4], fab[5], fab[8]+1))

            wx = np.hstack(wx)
            wy = np.hstack(wy)
            wz = np.hstack(wz)

            np.savetxt('wx_%i.txt' % (ilevel+1), zip(wx))
            np.savetxt('wy_%i.txt' % (ilevel+1), zip(wy))
            np.savetxt('wz_%i.txt' % (ilevel+1), zip(wz))
