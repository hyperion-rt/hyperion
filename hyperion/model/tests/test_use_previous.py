
from itertools import product

import pytest
import numpy as np

from .. import Model
from ...grid.amr_grid import AMRGrid, Level, Grid
from ...util.functions import random_id
from .test_helpers import get_test_dust


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize(('grid_type', 'copy'), list(product(['car', 'sph', 'cyl', 'amr', 'oct'], [(True, True), (False, True), (False, False)])))
def test_use_quantities(tmpdir, grid_type, copy):

    input_file_1 = tmpdir.join(random_id()).strpath
    output_file_1 = tmpdir.join(random_id()).strpath
    input_file_2 = tmpdir.join(random_id()).strpath
    output_file_2 = tmpdir.join(random_id()).strpath

    # Get a dust object
    dust = get_test_dust()
    dust_file = tmpdir.join(random_id()).strpath
    dust.write(dust_file)

    # Set up the initial model
    m = Model()

    # Set up initial geometry
    if grid_type == 'car':
        m.set_cartesian_grid([-1., 1.], [-2., 2.], [-3., 3.])
    elif grid_type == 'cyl':
        m.set_cylindrical_polar_grid([0., 1.], [-1., 1.], [0., 2. * np.pi])
    elif grid_type == 'sph':
        m.set_spherical_polar_grid([0., 1.], [0., np.pi], [0., 2. * np.pi])
    elif grid_type == 'amr':
        amr = AMRGrid()
        level = amr.add_level()
        grid = level.add_grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 8, 8
        grid.quantities['density'] = np.ones((8, 8, 8))
        m.set_amr_grid(amr)
    else:
        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        m.set_octree_grid(0., 0., 0., 10., 10., 10., np.array(refined).astype(bool))

    # Set up initial density
    if grid_type in ['car', 'cyl', 'sph']:
        m.add_density_grid(np.array([[[1.]]]), dust_file)
    elif grid_type == 'amr':
        m.add_density_grid(amr['density'], dust_file)
    else:
        m.add_density_grid(np.ones(len(refined)), dust_file)

    # Set up source
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 1000.

    # Set number of photons
    m.set_n_photons(initial=1000, imaging=1000)

    # Write out and run the model
    m.write(input_file_1, overwrite=True, copy=copy[1])
    m.run(output_file_1, overwrite=True)

    # Set up a second model that uses the properties from the first
    m2 = Model()

    # Use the geometry from the initial model
    m2.use_geometry(output_file_1)

    # Use the density and specific_energy from the initial model
    m2.use_quantities(output_file_1, copy=copy[0])

    # Set up source again
    s = m2.add_point_source()
    s.luminosity = 1.
    s.temperature = 1000.

    # Set number of photons
    m2.set_n_photons(initial=1000, imaging=1000)

    # Write out and run to test that the file is coherent
    m2.write(input_file_2, overwrite=True, copy=copy[1])
    m2.run(output_file_2, overwrite=True)


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize(('sublimation', 'output_density', 'output_density_diff'),
                         list(product(['no', 'fast'], ['none', 'last'], ['none', 'last'])))
def test_use_quantities_sublimation_warning(tmpdir, caplog, sublimation, output_density, output_density_diff):

    # If dust sublimation is enabled but the final density was not written
    # out, re-using the final quantities gives the input density, which may
    # not be consistent with the final specific energy, so a warning should
    # be emitted.

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    dust = get_test_dust()
    if sublimation != 'no':
        dust.set_sublimation_temperature(sublimation, temperature=1.)

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), dust)
    m.conf.output.output_density = output_density
    m.conf.output.output_density_diff = output_density_diff

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 1000.

    m.set_n_photons(initial=1000, imaging=0)
    m.set_n_initial_iterations(2)

    m.write(input_file)
    m.run(output_file)

    caplog.clear()
    Model.read(output_file, only_initial=True)
    assert 'The final density was not written' not in caplog.text

    caplog.clear()
    Model.read(output_file, only_initial=False)
    if sublimation != 'no' and output_density == 'none':
        assert 'The final density was not written' in caplog.text
        assert ('density_diff' in caplog.text) is (output_density_diff == 'last')
    else:
        assert 'The final density was not written' not in caplog.text
