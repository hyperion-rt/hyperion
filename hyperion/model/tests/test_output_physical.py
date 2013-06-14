from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np

from ...util.functions import random_id
from .test_helpers import get_test_model_noimaging, get_test_dust


@pytest.mark.parametrize(('output'), ['density', 'density_diff',
                                      'n_photons', 'specific_energy'])
def test_output_grids_exist(tmpdir, output):

    # Get a dust object
    dust = get_test_dust()

    # Set up the model
    model = get_test_model_noimaging()
    model.add_density_grid(np.array([[[1.]]]), dust)
    model.conf.output.output_density = 'last' if output == 'density' else 'none'
    model.conf.output.output_density_diff = 'last' if output == 'density_diff' else 'none'
    model.conf.output.output_n_photons = 'last' if output == 'n_photons' else 'none'
    model.conf.output.output_specific_energy = 'last' if output == 'specific_energy' else 'none'
    model.write(tmpdir.join(random_id()).strpath)

    # Run the model
    model_out = model.run(tmpdir.join(random_id()).strpath)

    # Check that component is available in output
    assert output in model_out.get_available_components()
    model_out.get_physical_grid(output)

    # If component is specific_energy, check that temperature is also available
    if output == 'specific_energy':
        assert 'temperature' in model_out.get_available_components()
        model_out.get_physical_grid('temperature')


def test_output_grids_density(tmpdir):

    # Get a dust object
    dust = get_test_dust()

    # Set initial density
    density_in = np.array([[[5.]]])

    # Set up the model
    model = get_test_model_noimaging()
    model.add_density_grid(density_in, dust)
    model.conf.output.output_density = 'last'
    model.conf.output.output_density_diff = 'none'
    model.conf.output.output_n_photons = 'none'
    model.conf.output.output_specific_energy = 'none'
    model.write(tmpdir.join(random_id()).strpath)

    # Run the model
    model_out = model.run(tmpdir.join(random_id()).strpath)

    # Extract density
    density_out = model_out.get_physical_grid('density', dust_id=0)
    assert density_in == density_out

    # Extract density (without specifying dust_id)
    density_out = model_out.get_physical_grid('density')
    assert density_in == density_out[0]
