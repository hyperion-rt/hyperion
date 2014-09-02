from __future__ import print_function, division

import numpy as np


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

    # TODO: generalize this once yt supports a custom refinement factor
    if refine != 2:
        raise ValueError("load_amr_grid only supports refine=2")

    xmin = ymin = zmin = +np.inf
    xmax = ymax = zmax = -np.inf

    grid_data = []

    for ilevel, level in enumerate(levels):

        for grid in level.grids:

            grid_dict = {}

            grid_dict['left_edge'] = [grid.xmin, grid.ymin, grid.zmin]
            grid_dict['right_edge'] = [grid.xmax, grid.ymax, grid.zmax]
            grid_dict['dimensions'] = [grid.nx, grid.ny, grid.nz]
            grid_dict['level'] = ilevel

            for field in grid.quantities:
                print(field, type(grid.quantities[field][dust_id]))
                grid_dict[field] = grid.quantities[field][dust_id]

            grid_data.append(grid_dict)

            xmin = min(xmin, grid.xmin)
            xmax = max(xmax, grid.xmax)
            ymin = min(ymin, grid.ymin)
            ymax = max(ymax, grid.ymax)
            zmin = min(zmin, grid.zmin)
            zmax = max(zmax, grid.zmax)

    # Determine domain resolution

    dx = (grid.xmax - grid.xmin) / float(grid.nx)
    nx = int(round((xmax - xmin) / dx))

    dy = (grid.ymax - grid.ymin) / float(grid.ny)
    ny = int(round((ymax - ymin) / dy))

    dz = (grid.zmax - grid.zmin) / float(grid.nz)
    nz = int(round((zmax - zmin) / dz))

    domain_dimensions = np.array([nx, ny, nz])

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    from yt.mods import load_amr_grids

    spf = load_amr_grids(grid_data, domain_dimensions, bbox=bbox)

    return spf


def octree_grid_to_yt_stream(grid, dust_id=0):

    xmin = grid.x - grid.dx
    xmax = grid.x + grid.dx
    ymin = grid.y - grid.dy
    ymax = grid.y + grid.dy
    zmin = grid.z - grid.dz
    zmax = grid.z + grid.dz

    from yt.mods import load_octree

    quantities = {}
    for field in grid.quantities:
        quantities[field] = grid.quantities[field][dust_id]

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    spf = load_octree(octree_mask=grid.refined.astype(np.uint8),
                      data=quantities,
                      bbox=bbox,
                      partial_coverage=1)

    return spf


def cartesian_grid_to_yt_stream(grid, xmin, xmax, ymin, ymax, zmin, zmax, dust_id=0):

    # TODO: only works for regular grids, need to catch non-uniform cases here

    # Make data dict which should contain (array, unit) tuples
    data = {}
    for field in grid.quantities:
        data[field] = (grid.quantities[field][dust_id], '')

    # Load cartesian grid into yt
    from yt.mods import load_uniform_grid
    spf = load_uniform_grid(data=data,
                            domain_dimensions=np.array(grid.shape[::-1], dtype=np.int32),
                            bbox=np.array([(xmin, xmax), (ymin, ymax), (zmin, zmax)]))

    return spf
