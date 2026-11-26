import numpy as np
import pytest
from numpy.testing import assert_allclose

from .. import Model, AnalyticalYSOModel
from ..helpers import tau_to_radius
from ...dust import IsotropicDust
from ...util.constants import c, rsun, lsun, tsun, au

# The models in this file have density that is constant along each radial
# column, so that the optical depth to escape radially from radius r is
# rho * chi * (r_max - r), and the radius at which it equals tau is
# r_max - tau / (rho * chi).

R_MAX = 10.


def constant_chi_dust(chi):
    return IsotropicDust([3.e9, 3.e16], [0., 0.], [chi, chi])


def spherical_model(n_r=100, n_theta=1):
    m = Model()
    m.set_spherical_polar_grid(np.linspace(0., R_MAX, n_r + 1),
                               np.linspace(0., np.pi, n_theta + 1),
                               [0., 2. * np.pi])
    return m


@pytest.mark.parametrize('tau', [0.1, 1., 3.3, 9.9])
def test_constant_density_1d(tau):

    rho, chi = 2., 0.5

    m = spherical_model()
    m.add_density_grid(np.ones(m.grid.shape) * rho, constant_chi_dust(chi))

    r = tau_to_radius(m, tau, 1.)

    assert r.shape == (1, 1)
    assert_allclose(r, R_MAX - tau / (rho * chi))


def test_tau_larger_than_total():

    # The total radial optical depth is 10, so larger values give zero

    m = spherical_model()
    m.add_density_grid(np.ones(m.grid.shape) * 2., constant_chi_dust(0.5))

    assert_allclose(tau_to_radius(m, 10.5, 1.), 0.)
    assert_allclose(tau_to_radius(m, 1.e3, 1.), 0.)


def test_wavelength_dependent_opacity():

    # Opacity proportional to frequency, which is interpolated exactly in
    # log-log space

    nu = np.array([3.e9, 3.e16])
    dust = IsotropicDust(nu, [0., 0.], nu / 3.e14)

    m = spherical_model()
    m.add_density_grid(np.ones(m.grid.shape), dust)

    for wav in [0.5, 1., 2.]:
        chi = c / (wav * 1.e-4) / 3.e14
        assert_allclose(tau_to_radius(m, 1., wav), R_MAX - 1. / chi)


def test_theta_dependent_density_2d():

    rho = np.array([0.2, 0.5, 1.0, 2.0])

    m = spherical_model(n_theta=len(rho))
    density = np.ones(m.grid.shape) * rho[np.newaxis, :, np.newaxis]
    m.add_density_grid(density, constant_chi_dust(1.))

    r = tau_to_radius(m, 1., 1.)

    assert r.shape == (1, len(rho))
    assert_allclose(r[0], R_MAX - 1. / rho)


def test_multiple_components():

    m = spherical_model()
    m.add_density_grid(np.ones(m.grid.shape) * 0.3, constant_chi_dust(1.))
    m.add_density_grid(np.ones(m.grid.shape) * 0.4, constant_chi_dust(0.5))

    # Total rho * chi is 0.3 + 0.2 = 0.5
    assert_allclose(tau_to_radius(m, 2., 1.), R_MAX - 2. / 0.5)


def test_not_model():
    with pytest.raises(TypeError) as exc:
        tau_to_radius(object(), 1., 1.)
    assert exc.value.args[0] == "model should be a Model instance"


def test_no_grid():
    with pytest.raises(Exception) as exc:
        tau_to_radius(Model(), 1., 1.)
    assert exc.value.args[0] == "Coordinate grid has not been defined"


def test_not_spherical():
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    with pytest.raises(TypeError) as exc:
        tau_to_radius(m, 1., 1.)
    assert exc.value.args[0] == "This method can only be called for spherical polar grids"

ENV_RHO_0 = 1.e-14
ENV_R_0 = au
ENV_R_MIN = au
ENV_R_MAX = 1000. * au
ENV_CHI = 100.


def power_law_envelope_model():

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun

    e = m.add_power_law_envelope()
    e.rmin = ENV_R_MIN
    e.rmax = ENV_R_MAX
    e.r_0 = ENV_R_0
    e.rho_0 = ENV_RHO_0
    e.power = -2.
    e.dust = IsotropicDust([3.e9, 3.e16], [0., 0.], [ENV_CHI, ENV_CHI])

    # to_model() copies the run configuration, which requires photon numbers
    m.set_n_photons(initial=100, imaging=100)

    return m


@pytest.mark.parametrize('tau', [0.01, 0.1, 1.])
def test_analytical_yso_auto_grid(tau):

    # For rho = rho_0 * (r / r_0) ** -2, the optical depth to escape radially
    # from r is rho_0 * r_0 ** 2 * chi * (1 / r - 1 / r_max)

    m = power_law_envelope_model()
    m.set_spherical_polar_grid_auto(1000, 1, 1)

    r = tau_to_radius(m, tau, 1.)

    expected = 1. / (tau / (ENV_RHO_0 * ENV_R_0 ** 2 * ENV_CHI) + 1. / ENV_R_MAX)

    assert r.shape == (1, 1)
    assert_allclose(r, expected, rtol=1.e-3)


def test_analytical_yso_explicit_grid():

    m = power_law_envelope_model()
    m.set_spherical_polar_grid(np.hstack([0., np.logspace(np.log10(ENV_R_MIN), np.log10(ENV_R_MAX), 1001)]),
                               [0., np.pi], [0., 2. * np.pi])

    r = tau_to_radius(m, 0.1, 1.)

    expected = 1. / (0.1 / (ENV_RHO_0 * ENV_R_0 ** 2 * ENV_CHI) + 1. / ENV_R_MAX)

    assert_allclose(r, expected, rtol=1.e-3)


def test_analytical_yso_auto_grid_cylindrical():
    m = power_law_envelope_model()
    m.set_cylindrical_polar_grid_auto(10, 10, 1)
    with pytest.raises(TypeError) as exc:
        tau_to_radius(m, 1., 1.)
    assert exc.value.args[0] == "This method can only be called for spherical polar grids"
