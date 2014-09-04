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

            grid_dict['left_edge'] = [grid.zmin, grid.ymin, grid.xmin]
            grid_dict['right_edge'] = [grid.zmax, grid.ymax, grid.xmax]
            grid_dict['dimensions'] = [grid.nz, grid.ny, grid.nx]
            grid_dict['level'] = ilevel

            for field in grid.quantities:
                grid_dict[('gas', field)] = grid.quantities[field][dust_id]

            grid_data.append(grid_dict)

            xmin = min(xmin, grid.xmin)
            xmax = max(xmax, grid.xmax)
            ymin = min(ymin, grid.ymin)
            ymax = max(ymax, grid.ymax)
            zmin = min(zmin, grid.zmin)
            zmax = max(zmax, grid.zmax)

    # Determine domain resolution

    grid0 = levels[0].grids[0]

    dx = (grid0.xmax - grid0.xmin) / float(grid0.nx)
    nx = int(round((xmax - xmin) / dx))

    dy = (grid0.ymax - grid0.ymin) / float(grid0.ny)
    ny = int(round((ymax - ymin) / dy))

    dz = (grid0.zmax - grid0.zmin) / float(grid0.nz)
    nz = int(round((zmax - zmin) / dz))

    domain_dimensions = np.array([nz, ny, nx])

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    from yt.mods import load_amr_grids

    spf = load_amr_grids(grid_data, domain_dimensions, bbox=bbox)

    return spf


def find_order(refined):
    """
    Find the index array to use to sort the ``refined`` and ``density`` arrays
    to swap the xyz <-> zyx order.
    """

    order = np.zeros(refined.shape)

    if not refined[0]:
        return [0]

    def find_nested(i):
        cells = [i]
        for cell in range(8):
            i += 1
            if refined[i]:
                parent = i
                i, sub_cells = find_nested(i)
                cells.append(sub_cells)
            else:
                cells.append(i)
        cells = [cells[j] for j in [0,1,5,3,7,2,6,4,8]]
        return i, np.hstack(cells)

    return find_nested(0)[1]


def octree_grid_to_yt_stream(grid, dust_id=0):

    order = find_order(grid.refined)
    refined = grid.refined[order]

    xmin = grid.x - grid.dx
    xmax = grid.x + grid.dx
    ymin = grid.y - grid.dy
    ymax = grid.y + grid.dy
    zmin = grid.z - grid.dz
    zmax = grid.z + grid.dz

    from yt.mods import load_octree

    quantities = {}
    for field in grid.quantities:
        quantities[('gas', field)] = np.atleast_2d(grid.quantities[field][dust_id][order][~refined]).transpose()

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    octree_mask = refined.astype(np.uint8) * 8

    spf = load_octree(octree_mask=octree_mask,
                      data=quantities,
                      bbox=bbox,
                      over_refine_factor=0,
                      partial_coverage=0)

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
                            domain_dimensions=np.array(grid.shape, dtype=np.int32),
                            bbox=np.array([(xmin, xmax), (ymin, ymax), (zmin, zmax)]))

    return spf
