import numpy as np

from ..amr_grid import AMRGrid, zero_density


def _single_grid_amr():
    amr = AMRGrid()
    level = amr.add_level()
    grid = level.add_grid()
    grid.xmin, grid.xmax = -1., 1.
    grid.ymin, grid.ymax = -1., 1.
    grid.zmin, grid.zmax = -1., 1.
    grid.nx, grid.ny, grid.nz = 4, 4, 4
    grid.quantities['density'] = np.ones((4, 4, 4))
    return amr


def test_zero_density_default_keeps_density():
    # Regression: zmin used to default to +inf, so the reset mask covered
    # every cell and the whole grid was zeroed. With no limits given, nothing
    # should be zeroed.
    amr = _single_grid_amr()
    out = zero_density(amr)
    assert out is amr                                  # returns the grid passed in
    assert np.all(out.levels[0].grids[0].quantities['density'] == 1.)


def test_zero_density_zeroes_outside_box():
    # Only cells with z > 0 should be zeroed.
    amr = _single_grid_amr()
    zero_density(amr, zmin=-np.inf, zmax=0.)
    dens = amr.levels[0].grids[0].quantities['density']
    # cell z-centres are at -0.75, -0.25, 0.25, 0.75 (first axis is z)
    assert np.all(dens[2:] == 0.)                      # z > 0 zeroed
    assert np.all(dens[:2] == 1.)                      # z < 0 kept
