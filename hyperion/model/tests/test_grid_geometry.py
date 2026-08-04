
import os
import tempfile

import numpy as np
import pytest

from .. import Model
from ...grid import AMRGrid
from ...util.functions import random_id
from .test_helpers import get_test_dust


@pytest.mark.requires_hyperion_binaries
def test_amr_source_outside_grid():
    # Regression test: for AMR grids, a photon emitted at a position outside
    # all level 1 grids used to index grids(-1) (undefined behavior, in
    # practice a segmentation fault). find_cell now returns invalid_cell, so
    # the run stops with a controlled, informative error instead.

    amr = AMRGrid()
    level = amr.add_level()
    grid = level.add_grid()
    grid.xmin, grid.xmax = -1., 1.
    grid.ymin, grid.ymax = -1., 1.
    grid.zmin, grid.zmax = -1., 1.
    grid.nx = grid.ny = grid.nz = 4
    grid.quantities['density'] = np.ones((4, 4, 4)) * 1.e-30

    m = Model()
    m.set_grid(amr)
    m.add_density_grid(amr['density'], get_test_dust())

    s = m.add_point_source()
    s.position = (2., 0., 0.)      # outside the level 1 grid
    s.luminosity = 1.
    s.temperature = 6000.

    m.set_n_initial_iterations(1)
    m.set_n_photons(initial=100, imaging=0)

    tmpdir = tempfile.mkdtemp()
    m.write(os.path.join(tmpdir, random_id()))
    logfile = os.path.join(tmpdir, 'log')

    with pytest.raises(SystemExit):
        m.run(os.path.join(tmpdir, 'out.rtout'), logfile=logfile)

    # The run should have failed through the controlled error path (with the
    # warning from find_cell and the explanatory error from emit), not
    # through a crash.
    log = open(logfile).read()
    assert 'photon not in any level 1 grid' in log
    assert 'not emitted inside a cell' in log


@pytest.mark.requires_hyperion_binaries
def test_theta_wall_viewing_angle():
    # Regression test for the theta cone wall "riding" condition: peeloff
    # rays whose direction has exactly the same polar angle as a theta wall
    # (e.g. an observer at 60 degrees with a wall at 60 degrees) but a
    # different azimuth used to be frozen on the wall when crossing it, and
    # were then killed, so the SED at that exact viewing angle collapsed to
    # zero. The source position and observer azimuth here are chosen so that
    # every direct peeloff ray at 60 degrees crosses the 60 degree cone.

    m = Model()
    r = np.array([0., 0.5, 1.0])
    t = np.array([0., np.radians(60.), np.pi])   # theta wall exactly at 60 deg
    p = np.array([0., 2. * np.pi])
    m.set_spherical_polar_grid(r, t, p)
    m.add_density_grid(np.ones((1, 2, 2)) * 1.e-8, get_test_dust())

    s = m.add_point_source()
    s.position = (0.3, 0., 0.05)    # equator side of the cone, azimuth 0
    s.luminosity = 1.
    s.temperature = 6000.

    sed = m.add_peeled_images(sed=True, image=False)
    # opposing azimuth ensures the peeloff rays cross the cone
    sed.set_viewing_angles([59., 60., 61.], [170., 170., 170.])
    sed.set_wavelength_range(10, 0.1, 100.)

    m.set_n_initial_iterations(0)
    # check the cell consistency at every step so that mis-binned photons
    # cannot go unnoticed
    m.set_propagation_check_frequency(1.)
    m.set_n_photons(imaging=5000)

    tmpdir = tempfile.mkdtemp()
    m.write(os.path.join(tmpdir, random_id()))
    mo = m.run(os.path.join(tmpdir, 'out.rtout'))

    wav, nufnu = mo.get_sed()
    nufnu = np.array(nufnu)
    tot = np.array([np.nansum(nufnu[i]) for i in range(3)])

    # all three flux totals should be positive and essentially identical:
    # the model is optically thin, so nothing special may happen at the
    # viewing angle that exactly matches the theta wall
    assert np.all(tot > 0.)
    np.testing.assert_allclose(tot[1], 0.5 * (tot[0] + tot[2]), rtol=0.1)
