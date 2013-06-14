from astropy.tests.helper import pytest
import numpy as np

from .. import Model
from .test_helpers import random_id, get_test_dust
from ...grid import AMRGrid


@pytest.mark.parametrize(('direction'), ['x', 'y', 'z'])
def test_amr_differing_widths(tmpdir, direction):

    # Widths of grids in same level are not the same

    dust = get_test_dust()

    amr = AMRGrid()

    level1 = amr.add_level()

    grid1 = level1.add_grid()
    grid1.nx = grid1.ny = grid1.nz = 4
    grid1.xmin = grid1.ymin = grid1.zmin = -10.
    grid1.xmax = grid1.ymax = grid1.zmax = +10.
    grid1.quantities['density'] = np.ones(grid1.shape) * 1.e-10

    grid2 = level1.add_grid()
    grid2.nx = grid2.ny = grid2.nz = 4
    grid2.xmin = grid2.ymin = grid2.zmin = -10.
    grid2.xmax = grid2.ymax = grid2.zmax = +10.
    grid2.quantities['density'] = np.ones(grid2.shape) * 1.e-10

    setattr(grid2, direction + 'min', -10.1)

    m = Model()
    m.set_grid(amr)
    m.add_density_grid(amr['density'], dust)
    m.set_n_photons(initial=1, imaging=0)
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)

    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert ('Grids 1 and 2 in level 1 have differing cell widths in the %s \n           direction ( 5.0000E+00 and  5.0250E+00 respectively)' % direction) in open(log_file).read()


@pytest.mark.parametrize(('direction'), ['x', 'y', 'z'])
def test_amr_misaligned_grids_same_level(tmpdir, direction):

    # Widths of grids in same level are not the same

    dust = get_test_dust()

    amr = AMRGrid()

    level1 = amr.add_level()

    grid1 = level1.add_grid()
    grid1.nx = grid1.ny = grid1.nz = 4
    grid1.xmin = grid1.ymin = grid1.zmin = -10.
    grid1.xmax = grid1.ymax = grid1.zmax = +10.
    grid1.quantities['density'] = np.ones(grid1.shape) * 1.e-10

    grid2 = level1.add_grid()
    grid2.nx = grid2.ny = grid2.nz = 4
    grid2.xmin = grid2.ymin = grid2.zmin = -10.
    grid2.xmax = grid2.ymax = grid2.zmax = +10.
    grid2.quantities['density'] = np.ones(grid2.shape) * 1.e-10

    setattr(grid2, direction + 'min', -10.1)
    setattr(grid2, direction + 'max', 9.9)

    m = Model()
    m.set_grid(amr)
    m.add_density_grid(amr['density'], dust)
    m.set_n_photons(initial=1, imaging=0)
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)

    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert ('Grids 1 and 2 in level 1 have edges that are not separated by \n           an integer number of cells in the %s direction' % direction) in open(log_file).read()


@pytest.mark.parametrize(('direction'), ['x', 'y', 'z'])
def test_amr_non_integer_refinement(tmpdir, direction):

    # Widths of grids in same level are not the same

    dust = get_test_dust()

    amr = AMRGrid()

    level1 = amr.add_level()

    grid1 = level1.add_grid()
    grid1.nx = grid1.ny = grid1.nz = 4
    grid1.xmin = grid1.ymin = grid1.zmin = -10.
    grid1.xmax = grid1.ymax = grid1.zmax = +10.
    grid1.quantities['density'] = np.ones(grid1.shape) * 1.e-10

    level2 = amr.add_level()

    grid2 = level2.add_grid()
    grid2.nx = grid2.ny = grid2.nz = 4
    grid2.xmin = grid2.ymin = grid2.zmin = -5.
    grid2.xmax = grid2.ymax = grid2.zmax = +5.
    grid2.quantities['density'] = np.ones(grid2.shape) * 1.e-10

    setattr(grid2, direction + 'min', -6.)

    m = Model()
    m.set_grid(amr)
    m.add_density_grid(amr['density'], dust)
    m.set_n_photons(initial=1, imaging=0)
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)

    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert ('Refinement factor in the %s direction between level 1 and \n           level 2 is not an integer (1.818)' % direction) in open(log_file).read()


@pytest.mark.parametrize(('direction'), ['x', 'y', 'z'])
def test_amr_not_aligned_across_levels(tmpdir, direction):

    # Widths of grids in same level are not the same

    dust = get_test_dust()

    amr = AMRGrid()

    level1 = amr.add_level()

    grid1 = level1.add_grid()
    grid1.nx = grid1.ny = grid1.nz = 4
    grid1.xmin = grid1.ymin = grid1.zmin = -10.
    grid1.xmax = grid1.ymax = grid1.zmax = +10.
    grid1.quantities['density'] = np.ones(grid1.shape) * 1.e-10

    level2 = amr.add_level()

    grid2 = level2.add_grid()
    grid2.nx = grid2.ny = grid2.nz = 4
    grid2.xmin = grid2.ymin = grid2.zmin = -5.
    grid2.xmax = grid2.ymax = grid2.zmax = +5.
    grid2.quantities['density'] = np.ones(grid2.shape) * 1.e-10

    setattr(grid2, direction + 'min', -6.)
    setattr(grid2, direction + 'max', 4.)

    m = Model()
    m.set_grid(amr)
    m.add_density_grid(amr['density'], dust)
    m.set_n_photons(initial=1, imaging=0)
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)

    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert ('Grid 1 in level 2 is not aligned with cells in level 1 in the \n           %s direction' % direction) in open(log_file).read()
