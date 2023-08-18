from __future__ import print_function, division

import numpy as np
import six
from astropy import log as logger


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
                grid_dict[('gas', field)] = grid.quantities[field][dust_id].transpose()

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

    domain_dimensions = np.array([nx, ny, nz])

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    from yt import load_amr_grids

    spf = load_amr_grids(grid_data, domain_dimensions, bbox=bbox,
                         geometry=('cartesian', ('x', 'y', 'z')))

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
        cells = [cells[j] for j in [0, 1, 5, 3, 7, 2, 6, 4, 8]]
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

    from yt import load_octree

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
        data[field] = grid.quantities[field][dust_id].transpose(), ''

    # Load cartesian grid into yt
    from yt import load_uniform_grid
    spf = load_uniform_grid(data=data,
                            domain_dimensions=np.array(grid.shape[::-1], dtype=np.int32),
                            bbox=np.array([(xmin, xmax), (ymin, ymax), (zmin, zmax)]),
                            geometry=('cartesian', ('x', 'y', 'z')))

    return spf


# For the following function, we use the low-level API described here:
# http://yt-project.org/docs/dev/examining/low_level_inspection.html


def yt_dataset_to_amr_grid(ds, quantity_mapping={}):
    """
    Convert a yt dataset to a Hyperion AMRGrid object

    .. note:: This function requires yt 3.0 or later

    Parameters
    ----------

    ds : yt Dataset
        The yt dataset
    quantity_mapping : dict
        A dictionary mapping the name of the quantity to use in Hyperion (the
        key) to the name of the field to extract in yt (the value). An example
        is provided below.

    Notes
    -----

    The domain is always re-centered so that the position at
    ds.domain_center in yt becomes the origin in Hyperion.

    Examples
    --------

    Assuming that your dust opacities are defined per unit gas mass, and the
    simulation density is given in gas densities, converting is
    straightfoward (in this case we assume the density field is called
    ``('gas', 'density')``)::

        >>> from yt import load
        >>> from hyperion.yt_wrappers import yt_dataset_to_amr_grid
        >>> ds = load('DD0010/moving7_0010')
        >>> amr = yt_dataset_to_amr_grid(ds, quantity_mapping={'density':('gas', 'density')})

    However, you will need to take care if your dust opacities are defined
    in dust mass units. If the yt dataset does not contain dust densities,
    you can add a field yourself, for example::

        >>> from yt import load
        >>> from hyperion.yt_wrappers import yt_dataset_to_amr_grid
        >>> ds = load('DD0010/moving7_0010')
        >>> def _dust_density(field, data):
        ...     return data[('gas', 'density')].in_units('g/cm**3') * 0.01
        >>> ds.add_field(('gas', 'dust_density'), function=_dust_density, units='g/cm**3', sampling_type='cell')

        >>> amr = yt_dataset_to_amr_grid(ds, quantity_mapping={'density':('gas', 'dust_density')})
    """

    field_list = "\n    ".join([str(x) for x in ds.derived_field_list])

    if len(quantity_mapping) == 0:
        raise ValueError("quantity_mapping needs to specified with key:value "
                         "pairs where the key is the name to give the quantity "
                         "in Hyperion and value is the name of the field in the "
                         "yt dataset. Available quantities are: \n\n    {0}".format(field_list))

    for output_quantity, input_field in six.iteritems(quantity_mapping):
        if not isinstance(output_quantity, six.string_types):
            raise ValueError("quantity_mapping keys should be strings")
        if input_field not in ds.derived_field_list:
            raise ValueError("yt field {0} does not exist. Available fields "
                             "are: \n\n    {1}".format(input_field, field_list))

    z0, y0, x0 = ds.domain_center.in_units('cm').ndarray_view()
    dz, dy, dx = ds.domain_width.in_units('cm').ndarray_view()

    logger.info("Domain center: x={0}cm, y={1}cm, z={2}cm".format(x0, y0, z0))
    logger.info("Domain width: dx={0}cm, dy={1}cm, dz={2}cm".format(dx, dy, dz))

    # Get levels and limits of all the grids
    n_levels = ds.index.max_level + 1
    levels = ds.index.grid_levels
    zmin, ymin, xmin = ds.index.grid_left_edge.in_units('cm').ndarray_view().transpose()
    zmax, ymax, xmax = ds.index.grid_right_edge.in_units('cm').ndarray_view().transpose()

    logger.info("Re-centering simulation so that domain center is at (0, 0, 0)")
    xmin -= x0
    xmax -= x0
    ymin -= y0
    ymax -= y0
    zmin -= z0
    zmax -= z0

    # Loop over levels and add grids
    from .amr_grid import AMRGrid
    amr = AMRGrid()
    for ilevel in range(n_levels):

        # Add a new level
        level = amr.add_level()

        # Loop over yt grids that are at this level
        for igrid in np.nonzero(levels == ilevel)[0]:

            # Get yt grid
            yt_grid = ds.index.grids[igrid]

            # Add a new Hyperion grid
            grid = level.add_grid()

            grid.xmin, grid.xmax = xmin[igrid], xmax[igrid]
            grid.ymin, grid.ymax = ymin[igrid], ymax[igrid]
            grid.zmin, grid.zmax = zmin[igrid], zmax[igrid]

            grid.nz, grid.ny, grid.nx = yt_grid.shape

            for output_quantity, input_field in six.iteritems(quantity_mapping):
                grid.quantities[output_quantity] = yt_grid[input_field].in_units('g/cm**3').ndarray_view()

    return amr
