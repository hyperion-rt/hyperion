import numpy as np

import yt.frontends.stream.api as stream
from yt.utilities.io_handler import BaseIOHandler


class HyperionIOHandler(BaseIOHandler):

    def __init__(self, grids, dust_id):
        if np.array(grids[0]).ndim == 4:
            self.dust_id = dust_id
        else:
            self.dust_id = None
        self.grids = grids
        BaseIOHandler.__init__(self)

    def _read_data_set(self, grid, field):
        if self.dust_id is None:
            return np.array(self.grids[grid.id][field.lower()].transpose())
        else:
            return np.array(self.grids[grid.id][field.lower()][self.dust_id].transpose())

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]


class StreamFieldData(object):

    def __init__(self, fields):
        self.fields = fields

    @property
    def all_fields(self):
        return self.fields


def almost_equal(a, b):
    return a / b < 1. + 1.e-4 and b / a < 1. + 1.e-4


def amr_grid_to_yt_stream(levels, dust_id=0):

    # Try and guess the refinement ratio - if it is not constant, then
    # we can't use yt

    if len(levels) == 0 or len(levels[0].grids) == 0:
        raise Exception("Need at least one level with one grid to convert to a yt object")
    elif len(levels) == 1:
        refine = 2
    else:
        dx = []
        dy = []
        dz = []
        for ilevel, level in enumerate(levels):
            for igrid, grid in enumerate(level.grids):
                gdx = (grid.xmax - grid.xmin) / float(grid.nx)
                gdy = (grid.ymax - grid.ymin) / float(grid.ny)
                gdz = (grid.zmax - grid.zmin) / float(grid.nz)
                if igrid == 0:
                    dx.append(gdx)
                    dy.append(gdy)
                    dz.append(gdz)
                else:
                    if not almost_equal(dx[-1], gdx):
                        raise Exception("dx scale differs between grids in level %i (expected %g and got %g)" % (ilevel, dx[-1], gdx))
                    if not almost_equal(dy[-1], gdy):
                        raise Exception("dy scale differs between grids in level %i (expected %g and got %g)" % (ilevel, dy[-1], gdy))
                    if not almost_equal(dz[-1], gdz):
                        raise Exception("dz scale differs between grids in level %i (expected %g and got %g)" % (ilevel, dz[-1], gdz))
        dx = np.array(dx)
        dy = np.array(dy)
        dz = np.array(dz)
        refine_x = dx[:-1] / dx[1:]
        refine_y = dy[:-1] / dy[1:]
        refine_z = dz[:-1] / dz[1:]
        for i in range(len(levels) - 1):
            if abs(refine_x[i] - round(refine_x[i])) > 1.e-5:
                raise Exception("refinement ratio is not an integer (%g)" % refine_x[i])
            if abs(refine_y[i] - round(refine_y[i])) > 1.e-5:
                raise Exception("refinement ratio is not an integer (%g)" % refine_y[i])
            if abs(refine_z[i] - round(refine_z[i])) > 1.e-5:
                raise Exception("refinement ratio is not an integer (%g)" % refine_z[i])
        refine_x = np.round(refine_x).astype(int)
        refine_y = np.round(refine_y).astype(int)
        refine_z = np.round(refine_z).astype(int)
        if not np.all(np.hstack([refine_x, refine_y, refine_z]) == refine_x[0]):
            raise Exception("refinement ratio changes between levels and/or directions (x = %s, y = %s, z = %s)" % (str(refine_x), str(refine_y), str(refine_z)))
        refine = int(refine_x[0])

    xmin = ymin = zmin = +np.inf
    xmax = ymax = zmax = -np.inf

    left_edge = []
    right_edge = []
    dimensions = []
    level_ids = []

    n_grids = 0
    grids = []

    for ilevel, level in enumerate(levels):

        for grid in level.grids:

            left_edge.append([grid.xmin,
                              grid.ymin,
                              grid.zmin])

            right_edge.append([grid.xmax,
                               grid.ymax,
                               grid.zmax])

            dimensions.append([grid.nx,
                               grid.ny,
                               grid.nz])

            level_ids.append(ilevel)

            n_grids += 1

            xmin = min(xmin, grid.xmin)
            xmax = max(xmax, grid.xmax)
            ymin = min(ymin, grid.ymin)
            ymax = max(ymax, grid.ymax)
            zmin = min(zmin, grid.zmin)
            zmax = max(zmax, grid.zmax)

            grids.append(grid.quantities)

    left_edge = np.array(left_edge, dtype=np.float64)
    right_edge = np.array(right_edge, dtype=np.float64)
    dimensions = np.array(dimensions, dtype=np.int32)

    # The reshape is necessary due to a strange Numpy bug
    level_ids = np.array(level_ids, dtype=np.int32).reshape((len(level_ids), 1))

    parent_ids = None
    particle_count = np.zeros((n_grids, 1), dtype='int32')
    processor_ids = np.zeros(n_grids)

    grid = levels[0].grids[0]

    # Determine domain resolution

    dx = (grid.xmax - grid.xmin) / float(grid.nx)
    nx = int(round((xmax - xmin) / dx))

    dy = (grid.ymax - grid.ymin) / float(grid.ny)
    ny = int(round((ymax - ymin) / dy))

    dz = (grid.zmax - grid.zmin) / float(grid.nz)
    nz = int(round((zmax - zmin) / dz))

    # Determine fields

    fields = grid.quantities.keys()

    # Set up StreamHandler

    handler = stream.StreamHandler(
        left_edge[:],
        right_edge[:],
        dimensions[:],
        level_ids[:],
        parent_ids,
        particle_count[:],
        processor_ids[:],
        StreamFieldData(fields),
        HyperionIOHandler(grids, dust_id),
    )

    handler.name = 'hyperion'
    handler.domain_left_edge = np.array([xmin, ymin, zmin])
    handler.domain_right_edge = np.array([xmax, ymax, zmax])
    handler.refine_by = refine
    handler.dimensionality = 3
    handler.domain_dimensions = np.array([nx, ny, nz])
    handler.simulation_time = 0.0
    handler.cosmology_simulation = 0

    spf = stream.StreamStaticOutput(handler)
    spf.units["cm"] = 1.0
    spf.units["unitary"] = 1.0 / ((spf.domain_right_edge - spf.domain_left_edge).max())

    return spf
