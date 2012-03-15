import numpy as np

from hyperion.util.meshgrid import meshgrid_nd
from hyperion.util.functions import FreezableClass


def zero_density(grid, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=np.inf, zmax=np.inf):
    for ilevel,level in enumerate(grid.levels):
        for igrid,grid in enumerate(level.grids):
            wx = np.linspace(grid.xmin, grid.xmax, grid.nx + 1)
            wy = np.linspace(grid.ymin, grid.ymax, grid.ny + 1)
            wz = np.linspace(grid.zmin, grid.zmax, grid.nz + 1)
            x = 0.5 * (wx[:-1] + wx[1:])
            y = 0.5 * (wy[:-1] + wy[1:])
            z = 0.5 * (wz[:-1] + wz[1:])
            gx, gy, gz = meshgrid_nd(x, y, z)
            reset = (gx < xmin) | (gx > xmax) | (gy < ymin) | (gy > ymax) | (gz < zmin) | (gz > zmax)
            grid.data[reset] = 0.
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
            for igrid, grid in enumerate(level.grids):
                grid_name = "Grid %i" % (igrid + 1)
                if grid_name in g_level:
                    g_grid = g_level[grid_name]
                else:
                    g_grid = g_level.create_group("Grid %i" % (igrid + 1))
                if dust:

                    # First check that the dimensions are ok
                    for a in array:
                        if a.levels[ilevel].grids[igrid].data.shape != (grid.nz,grid.ny,grid.nx):
                            raise Exception("Grid dimensions inconsistent with physical array dimensions")

                    # Restack list of arrays as single array
                    shape = list(array[0].levels[ilevel].grids[igrid].data.shape)
                    shape.insert(0, len(array))
                    grid_array = np.vstack([a.levels[ilevel].grids[igrid].data for a in array]).reshape(*shape)

                else:

                    # First check that the dimensions are ok
                    if array.levels[ilevel].grids[igrid].data.shape != (grid.nz,grid.ny,grid.nx):
                        raise Exception("Grid dimensions inconsistent with physical array dimensions")

                    grid_array = array.levels[ilevel].grids[igrid].data

                g_grid.create_dataset(name, data=grid_array, compression=compression, dtype=physics_dtype)

    def write_geometry(self, group, overwrite=True, volumes=False, areas=False, widths=False, compression=True, geo_dtype=float, wall_dtype=float):
        group.attrs['geometry'] = self.geometry_id
        group.attrs['grid_type'] = 'amr'
        group.attrs['nlevels'] = len(self.levels)
        for ilevel, level in enumerate(self.levels):
            g_level = group.create_group("Level %i" % (ilevel + 1))
            g_level.attrs['ngrids'] = len(level.grids)
            for igrid, grid in enumerate(level.grids):
                g_grid = g_level.create_group("Grid %i" % (igrid + 1))
                g_grid.attrs['xmin'] = grid.xmin
                g_grid.attrs['xmax'] = grid.xmax
                g_grid.attrs['ymin'] = grid.ymin
                g_grid.attrs['ymax'] = grid.ymax
                g_grid.attrs['zmin'] = grid.zmin
                g_grid.attrs['zmax'] = grid.zmax
                g_grid.attrs['n1'] = grid.nx
                g_grid.attrs['n2'] = grid.ny
                g_grid.attrs['n3'] = grid.nz
