import os

import pytest
import numpy as np

from .. import Model, AnalyticalYSOModel, ModelOutput
from ..helpers import hseq_profile, run_with_vertical_hseq
from ...util.constants import G, m_h, k, msun, rsun, lsun, tsun, au
from ...util.integrate import integrate
from .test_helpers import get_test_dust

MU = 2.279


def isothermal_profile(w, z, temperature, mstar, mu=MU):
    factor = G * mstar * mu * m_h / (k * temperature)
    density = np.exp(-factor * (1. / w - 1. / np.sqrt(w ** 2 + z ** 2)))
    return density / integrate(z, density)


def scale_height(w, temperature, mstar, mu=MU):
    return np.sqrt(k * temperature * w ** 3 / (G * mstar * mu * m_h))


# Use an even number of points so that z=0 is not in the grid, as for cell
# centres in a grid symmetric about the mid-plane
@pytest.mark.parametrize('n_z', [400, 401])
def test_hseq_profile_isothermal(n_z):

    w = au
    temperature = 100.
    h = scale_height(w, temperature, msun)

    z = np.linspace(-6. * h, 6. * h, n_z)

    density = hseq_profile(w, z, np.repeat(temperature, n_z), msun)

    # Normalized
    np.testing.assert_allclose(integrate(z, density), 1., rtol=1.e-10)

    # Symmetric about the mid-plane
    np.testing.assert_allclose(density, density[::-1], rtol=1.e-10)

    # Matches the analytical isothermal solution
    expected = isothermal_profile(w, z, temperature, msun)
    np.testing.assert_allclose(density, expected, rtol=1.e-3, atol=1.e-6 * expected.max())

    # Approximately Gaussian for z << w
    gaussian = np.exp(-0.5 * (z / h) ** 2)
    gaussian /= integrate(z, gaussian)
    np.testing.assert_allclose(density, gaussian, rtol=0.1, atol=1.e-3 * gaussian.max())


def test_hseq_profile_temperature_dependence():

    # A hotter column should be more extended vertically
    w = au
    z = np.linspace(-0.5 * au, 0.5 * au, 400)

    cold = hseq_profile(w, z, np.repeat(30., z.size), msun)
    hot = hseq_profile(w, z, np.repeat(300., z.size), msun)

    def rms_height(density):
        return np.sqrt(integrate(z, z ** 2 * density))

    np.testing.assert_allclose(rms_height(hot) / rms_height(cold), np.sqrt(10.), rtol=0.05)


def hseq_model():

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun
    m.star.mass = msun

    d = m.add_flared_disk()
    d.mass = 1.e-3 * msun
    d.rmin = au
    d.r_0 = au
    d.rmax = 10. * au
    d.h_0 = 0.05 * au
    d.p = -1.
    d.beta = 1.25
    d.dust = get_test_dust()

    m.set_cylindrical_polar_grid_auto(10, 20, 1)

    m.set_n_photons(initial=10000, imaging=0)
    m.set_n_initial_iterations(3)

    return m


def test_hseq_invalid_model(tmpdir):
    m = Model()
    m.set_cylindrical_polar_grid([0., au], [-au, au], [0., 2. * np.pi])
    with pytest.raises(TypeError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Can only run hydrostatic equilibrium for AnalyticalYSOModel instances'


def test_hseq_invalid_grid_auto(tmpdir):
    m = hseq_model()
    m.set_spherical_polar_grid_auto(10, 10, 1)
    with pytest.raises(TypeError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Can only run hydrostatic equilibrium for models with cylindrical polar grids'


@pytest.mark.xfail(raises=KeyError, strict=True,
                   reason="grid type check only works for grids set with set_*_grid_auto")
def test_hseq_invalid_grid_explicit(tmpdir):
    m = hseq_model()
    m.set_spherical_polar_grid([0., au, 2. * au], [0., np.pi], [0., 2. * np.pi])
    with pytest.raises(TypeError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Can only run hydrostatic equilibrium for models with cylindrical polar grids'


@pytest.mark.xfail(raises=KeyError, strict=True,
                   reason="grid type check only works for grids set with set_*_grid_auto")
def test_hseq_valid_grid_explicit(tmpdir):
    # An explicit cylindrical polar grid should pass the grid check, so we
    # unset the stellar mass to stop at the next check
    m = hseq_model()
    m.star.mass = None
    m.set_cylindrical_polar_grid([0., au, 2. * au], [-au, au], [0., 2. * np.pi])
    with pytest.raises(ValueError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Stellar mass needs to be defined for calculation of hydrostatic equilibrium'


def test_hseq_no_mass(tmpdir):
    m = hseq_model()
    m.star.mass = None
    with pytest.raises(ValueError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Stellar mass needs to be defined for calculation of hydrostatic equilibrium'


def test_hseq_no_disks(tmpdir):
    m = hseq_model()
    m.disks = []
    with pytest.raises(ValueError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m)
    assert exc.value.args[0] == 'Can only run hydrostatic equilibrium for models with disks'


def test_hseq_invalid_logfile(tmpdir):
    m = hseq_model()
    with pytest.raises(TypeError) as exc:
        run_with_vertical_hseq(tmpdir.join('model').strpath, m, logfile='test.log')
    assert exc.value.args[0] == 'logfile should be True or False'


@pytest.mark.requires_hyperion_binaries
def test_hseq_run(tmpdir):

    m = hseq_model()

    prefix = tmpdir.join('model').strpath

    run_with_vertical_hseq(prefix, m, n_iter=2, logfile=True)

    for iteration in range(3):
        for extension in ('rtin', 'rtout', 'log'):
            assert os.path.exists('{0:s}_{1:05d}.{2:s}'.format(prefix, iteration, extension))

    previous = None

    for iteration in range(3):

        g = ModelOutput('{0:s}_{1:05d}.rtout'.format(prefix, iteration)).get_quantities()

        density = g['density'][0].array

        assert np.all(np.isfinite(density))
        assert np.all(density >= 0.)

        # The surface density integrated over z should be conserved
        sigma = np.array([[integrate(g.z, density[j, :, i])
                           for i in range(len(g.w))]
                          for j in range(len(g.p))])

        if previous is None:
            assert np.any(sigma > 0.)
        else:
            np.testing.assert_allclose(sigma, previous[0], rtol=1.e-10)
            assert not np.allclose(density, previous[1], rtol=1.e-3, atol=0.)

        previous = sigma, density
