import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass


def zero_density(grid, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=np.inf, zmax=np.inf):
    for ilevel,level in enumerate(grid.levels):
        for ifab,fab in enumerate(level.fabs):
            wx = np.linspace(fab.xmin, fab.xmax, fab.nx + 1)
            wy = np.linspace(fab.ymin, fab.ymax, fab.ny + 1)
            wz = np.linspace(fab.zmin, fab.zmax, fab.nz + 1)
            x = 0.5 * (wx[:-1] + wx[1:])
            y = 0.5 * (wy[:-1] + wy[1:])
            z = 0.5 * (wz[:-1] + wz[1:])
            gx, gy, gz = meshgrid_nd(x, y, z)
            reset = (gx < xmin) | (gx > xmax) | (gy < ymin) | (gy > ymax) | (gz < zmin) | (gz > zmax)
            fab.data[reset] = 0.
    return grid


class AMRGrid(FreezableClass):

    def __init__(self, grid):

        self.levels = grid.levels

        self.geometry_id = None

        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'shape':
            return (1, 1, self.ncells)

    def write_physical_array(self, group, array, name, dust=False, compression=True, physics_dtype=float):

        for ilevel, level in enumerate(self.levels):
            level_name = "Level %i" % (ilevel + 1)
            if level_name in group:
                g_level = group[level_name]
            else:
                g_level = group.create_group("Level %i" % (ilevel + 1))
            for ifab, fab in enumerate(level.fabs):
                fab_name = "Fab %i" % (ifab + 1)
                if fab_name in g_level:
                    g_fab = g_level[fab_name]
                else:
                    g_fab = g_level.create_group("Fab %i" % (ifab + 1))
                if dust:

                    # First check that the dimensions are ok
                    for a in array:
                        if a.levels[ilevel].fabs[ifab].data.shape != (fab.nz,fab.ny,fab.nx):
                            raise Exception("Fab dimensions inconsistent with physical array dimensions")

                    # Restack list of arrays as single array
                    shape = list(array[0].levels[ilevel].fabs[ifab].data.shape)
                    shape.insert(0, len(array))
                    fab_array = np.vstack([a.levels[ilevel].fabs[ifab].data for a in array]).reshape(*shape)

                else:

                    # First check that the dimensions are ok
                    if array.levels[ilevel].fabs[ifab].data.shape != (fab.nz,fab.ny,fab.nx):
                        raise Exception("Fab dimensions inconsistent with physical array dimensions")

                    fab_array = array.levels[ilevel].fabs[ifab].data

                g_fab.create_dataset(name, data=fab_array, compression=compression, dtype=physics_dtype)

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):
        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'amr'
        group.attrs['nlevels'] = len(self.levels)
        for ilevel, level in enumerate(self.levels):
            g_level = group.create_group("Level %i" % (ilevel + 1))
            g_level.attrs['nfabs'] = len(level.fabs)
            for ifab, fab in enumerate(level.fabs):
                g_fab = g_level.create_group("Fab %i" % (ifab + 1))
                g_fab.attrs['xmin'] = fab.xmin
                g_fab.attrs['xmax'] = fab.xmax
                g_fab.attrs['ymin'] = fab.ymin
                g_fab.attrs['ymax'] = fab.ymax
                g_fab.attrs['zmin'] = fab.zmin
                g_fab.attrs['zmax'] = fab.zmax
                g_fab.attrs['n1'] = fab.nx
                g_fab.attrs['n2'] = fab.ny
                g_fab.attrs['n3'] = fab.nz
