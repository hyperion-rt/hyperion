from astropy.tests.helper import pytest
import numpy as np

from ..model import Model
from .test_helpers import get_realistic_test_dust, random_id

# The following tests ensure that the correct temperatures are returned
# for the MRW, for single and for multiple dust populations.

D_REF = np.logspace(-5., 12., 18)

T_REF = [24.75280,
         24.66414,
         24.52175,
         21.97109,
         15.53059,
         10.76363,
         7.810127,
         6.672520,
         6.902798,
         16.64318,
         53.08394,
         162.0158,
         438.9026,
         1013.141,
         2156.520,
         4642.825,
         9948.065,
         21211.08]


@pytest.mark.parametrize(('density_ref', 'temperature_ref'), zip(D_REF, T_REF))
def test_single_temperature(tmpdir, density_ref, temperature_ref):

    dust = get_realistic_test_dust()

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    source = m.add_point_source()
    source.luminosity = 1.
    source.temperature = 6000.

    m.add_density_grid(np.ones(m.grid.shape) * density_ref, dust)

    m.set_n_photons(initial=1000, imaging=0)

    m.set_n_initial_iterations(30)
    m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

    m.set_mrw(True, gamma=2)

    m.set_max_interactions(1000000000)

    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)

    grid = mo.get_quantities()
    temperature = grid['temperature'][0].array[0, 0, 0]

    assert temperature_ref / temperature < 1.1 and temperature / temperature_ref < 1.1


@pytest.mark.parametrize(('density_ref', 'temperature_ref'), zip(D_REF, T_REF))
def test_multi_temperature(tmpdir, density_ref, temperature_ref):

    dust = get_realistic_test_dust()

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    source = m.add_point_source()
    source.luminosity = 1.
    source.temperature = 6000.

    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.1, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.2, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.3, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.4, dust)

    m.set_n_photons(initial=1000, imaging=0)

    m.set_n_initial_iterations(30)
    m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

    m.set_mrw(True, gamma=2)

    m.set_max_interactions(1000000000)

    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)

    grid = mo.get_quantities()

    for i in range(4):
        temperature = grid['temperature'][i].array[0, 0, 0]
        assert temperature_ref / temperature < 1.1 and temperature / temperature_ref < 1.1
