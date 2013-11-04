from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np

from .. import AnalyticalYSOModel
from ...util.constants import msun, rsun, lsun, tsun, au, yr
from ...util.functions import random_id
from .test_helpers import get_test_dust
from ...util.convenience import OptThinRadius


def basic_analytical_model():

    dust = get_test_dust()

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun

    d = m.add_flared_disk()
    d.mass = 1.e-2 * msun
    d.rmin = 0.1 * au
    d.r_0 = au
    d.rmax = 300. * au
    d.h_0 = 0.01 * au
    d.p = -1.
    d.beta = 1.25
    d.dust = dust

    return m


def test_analytical_yso_full(tmpdir):

    m = basic_analytical_model()

    m.set_spherical_polar_grid_auto(10, 10, 10)

    m.set_n_photons(initial=100, imaging=100)

    m.write(tmpdir.join(random_id()).strpath)


def test_analytical_yso_nogrid_invalid(tmpdir):

    m = AnalyticalYSOModel()
    with pytest.raises(Exception) as e:
        m.write(tmpdir.join(random_id()).strpath)
    assert e.value.args[0] == 'The coordinate grid needs to be defined'


def test_analytical_yso_nostar_invalid():

    m = AnalyticalYSOModel()
    m.set_spherical_polar_grid_auto(1, 1, 1)
    with pytest.raises(Exception) as e:
        m.to_model()
    assert e.value.args[0] == 'The central source radius need to be defined before the grid can be set up'


def test_analytical_yso_optthinradius(tmpdir):

    dust = get_test_dust()

    m = basic_analytical_model()

    e = m.add_power_law_envelope()
    e.mass = 1.e-2 * msun
    e.rmin = OptThinRadius(1000.)
    e.r_0 = au
    e.rmax = OptThinRadius(10.)
    e.power = -2.
    e.dust = dust

    m.set_spherical_polar_grid_auto(10, 10, 10)

    m.set_n_photons(initial=100, imaging=100)

    m.write(tmpdir.join(random_id()).strpath)


def test_analytical_yso_optthinradius_check_not_frozen(tmpdir):
    """
    This test used to ensure that stellar parameters could no longer be
    changed after the grid was set, but now that the rmin/rmax properties of
    the density components are dynamic, this is no longer needed, so we
    instead check that the properties *can* be updated.
    """

    dust = get_test_dust()

    m = basic_analytical_model()

    e = m.add_power_law_envelope()
    e.mass = 1.e-2 * msun
    e.rmin = OptThinRadius(1000.)
    e.r_0 = au
    e.rmax = OptThinRadius(10.)
    e.power = -2.
    e.dust = dust

    m.set_spherical_polar_grid_auto(10, 10, 10)

    e.mass = 1.
    m.star.radius = 1.

    m.set_n_photons(initial=100, imaging=100)

    m.write(tmpdir.join(random_id()).strpath)


def test_analytical_yso_optthinradius_check_dynamic_evaluation():

    dust = get_test_dust()

    m = basic_analytical_model()

    e = m.add_power_law_envelope()
    e.mass = 1.e-2 * msun
    e.rmin = OptThinRadius(1000.)
    e.r_0 = au
    e.rmax = OptThinRadius(10.)
    e.power = -2.
    e.dust = dust

    assert np.isscalar(e.rmin)
    assert np.isscalar(e.rmax)


def test_analytical_yso_add_density():
    m = AnalyticalYSOModel()
    with pytest.raises(NotImplementedError) as exc:
        m.add_density_grid()
    assert exc.value.args[0] == 'add_density_grid cannot be used for AnalyticalYSOModel'


def test_analytical_yso_use_quantities_invalid(tmpdir):

    output_file = tmpdir.join(random_id()).strpath

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(tmpdir.join(random_id()).strpath)
    m.run(output_file)

    m2 = basic_analytical_model()
    m2.set_spherical_polar_grid_auto(10, 10, 10)
    m2.set_n_photons(initial=100, imaging=100)
    with pytest.raises(NotImplementedError) as exc:
        m2.use_quantities(output_file)
    assert exc.value.args[0] == "use_quantities cannot be used for AnalyticalYSOModel"


def test_analytical_yso_use_geometry_invalid(tmpdir):

    output_file = tmpdir.join(random_id()).strpath

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(tmpdir.join(random_id()).strpath)
    m.run(output_file)

    m2 = basic_analytical_model()
    with pytest.raises(NotImplementedError) as exc:
        m2.use_geometry(output_file)
    assert exc.value.args[0] == "use_geometry cannot be used for AnalyticalYSOModel"


def test_ambient_medium(tmpdir):

    dust = get_test_dust()

    m = AnalyticalYSOModel()

    m.star.radius = 1.
    m.star.temperature = 1000.
    m.star.luminosity = 1.

    d = m.add_flared_disk()
    d.mass = 1.
    d.rmin = 0.1
    d.rmax = 10.
    d.p = -1
    d.beta = 1.25
    d.h_0 = 0.1
    d.r_0 = 2.
    d.dust = dust

    a = m.add_ambient_medium()
    a.rmin = 0.1
    a.rmax = 10.
    a.rho = 1.
    a.dust = dust

    m.set_spherical_polar_grid_auto(399, 199, 1)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)


@pytest.mark.parametrize(('grid_type'), ['cylindrical', 'spherical'])
def test_rmin_zero(grid_type):

    m = AnalyticalYSOModel()

    m.star.radius = 1.
    m.star.temperature = 1000.
    m.star.luminosity = 1.

    d = m.add_flared_disk()
    d.rmin = 0.
    d.rmax = 10.
    d.r_0 = 10.
    d.h_0 = 1.
    d.p = -1
    d.beta = 1.25
    d.dust = get_test_dust()

    if grid_type == 'cylindrical':

        m.set_cylindrical_polar_grid_auto(100, 20, 3)

        with pytest.raises(ValueError) as exc:
            m.to_model()
        assert exc.value.args[0] == "R_min is 0, so cannot set up the grid cell walls automatically. Use set_cylindrical_polar_grid() instead to specify the cell wall positions."

    else:

        m.set_spherical_polar_grid_auto(100, 20, 3)

        with pytest.raises(ValueError) as exc:
            m.to_model()
        assert exc.value.args[0] == "R_min is 0, so cannot set up the grid cell walls automatically. Use set_spherical_polar_grid() instead to specify the cell wall positions."


def test_complete_spherical(tmpdir):

    m = AnalyticalYSOModel()

    m.star.radius = 0.01
    m.star.temperature = 1000.
    m.star.luminosity = 1.
    m.star.mass = 1.

    d = m.add_alpha_disk()
    d.mass = 0.001
    d.rmin = 0.1
    d.rmax = 10.
    d.r_0 = 10.
    d.h_0 = 1.
    d.p = -1
    d.beta = 1.25
    d.mdot = 1.
    d.dust = get_test_dust()

    d = m.add_flared_disk()
    d.mass = 0.001
    d.rmin = 0.1
    d.rmax = 10.
    d.r_0 = 10.
    d.h_0 = 1.
    d.p = -1
    d.beta = 1.25
    d.dust = get_test_dust()

    e = m.add_ulrich_envelope()
    e.rmin = 0.1
    e.rmax = 10.
    e.rho_0 = 1.
    e.rc = 2.
    e.dust = get_test_dust()

    c1 = e.add_bipolar_cavity()
    c1.dust = get_test_dust()
    c1.theta_0 = 10.
    c1.power = 1.2
    c1.r_0 = 3.
    c1.rho_0 = 2.3

    p = m.add_power_law_envelope()
    p.rmin = 0.2
    p.rmax = 11.
    p.mass = 10.
    p.power = -1.
    p.r_0 = 2.
    p.dust = get_test_dust()

    c2 = p.add_bipolar_cavity()
    c2.dust = get_test_dust()
    c2.theta_0 = 20
    c2.power = 1.3
    c2.r_0 = 4.
    c2.rho_0 = 9.9

    a = m.add_ambient_medium()
    a.rmin = 0.1
    a.rmax = 10.
    a.rho = 1.
    a.dust = get_test_dust()

    m.set_spherical_polar_grid_auto(399, 199, 1)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)


def test_complete_cylindrical(tmpdir):

    m = AnalyticalYSOModel()

    m.star.radius = 0.01
    m.star.temperature = 1000.
    m.star.luminosity = 1.
    m.star.mass = 1.

    d = m.add_alpha_disk()
    d.mass = 0.001
    d.rmin = 0.1
    d.rmax = 10.
    d.r_0 = 10.
    d.h_0 = 1.
    d.p = -1
    d.beta = 1.25
    d.mdot = 1.
    d.dust = get_test_dust()

    d = m.add_flared_disk()
    d.mass = 0.001
    d.rmin = 0.1
    d.rmax = 10.
    d.r_0 = 10.
    d.h_0 = 1.
    d.p = -1
    d.beta = 1.25
    d.dust = get_test_dust()

    m.set_cylindrical_polar_grid_auto(399, 199, 1)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)


def test_complete_spherical_optthin(tmpdir):

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun
    m.star.mass = msun

    d = m.add_alpha_disk()
    d.mass = 0.001 * msun
    d.rmin = OptThinRadius(1600.)
    d.rmax = OptThinRadius(30.)
    d.r_0 = 10. * au
    d.h_0 = 1. * au
    d.p = -1
    d.beta = 1.25
    d.mdot = 1.e-10 * msun / yr
    d.dust = get_test_dust()

    d = m.add_flared_disk()
    d.mass = 0.001 * msun
    d.rmin = OptThinRadius(1500.)
    d.rmax = OptThinRadius(20.)
    d.r_0 = 10. * au
    d.h_0 = 1. * au
    d.p = -1
    d.beta = 1.25
    d.dust = get_test_dust()

    e = m.add_ulrich_envelope()
    e.rmin = OptThinRadius(1400.)
    e.rmax = OptThinRadius(30.)
    e.rho_0 = 1.e-20
    e.rc = 2. * au
    e.dust = get_test_dust()

    c1 = e.add_bipolar_cavity()
    c1.dust = get_test_dust()
    c1.theta_0 = 10.
    c1.power = 1.2
    c1.r_0 = 3. * au
    c1.rho_0 = 2.3e-23

    p = m.add_power_law_envelope()
    p.rmin = OptThinRadius(500)
    p.rmax = OptThinRadius(20.)
    p.mass = 10. * msun
    p.power = -1.
    p.r_0 = 2. * au
    p.dust = get_test_dust()

    c2 = p.add_bipolar_cavity()
    c2.dust = get_test_dust()
    c2.theta_0 = 20
    c2.power = 1.3
    c2.r_0 = 4. * au
    c2.rho_0 = 9.9e-22

    a = m.add_ambient_medium()
    a.rmin = OptThinRadius(1600.)
    a.rmax = OptThinRadius(10.)
    a.rho = 1.e-18
    a.dust = get_test_dust()

    m.set_spherical_polar_grid_auto(399, 199, 1)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)


def test_complete_cylindrical_optthin(tmpdir):

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun
    m.star.mass = msun

    d = m.add_alpha_disk()
    d.mass = 0.001 * msun
    d.rmin = OptThinRadius(1600.)
    d.rmax = OptThinRadius(30.)
    d.r_0 = 10. * au
    d.h_0 = 1. * au
    d.p = -1
    d.beta = 1.25
    d.mdot = 1.e-10 * msun / yr
    d.dust = get_test_dust()

    d = m.add_flared_disk()
    d.mass = 0.001 * au
    d.rmin = OptThinRadius(1500.)
    d.rmax = OptThinRadius(20.)
    d.r_0 = 10. * au
    d.h_0 = 1. * au
    d.p = -1
    d.beta = 1.25
    d.dust = get_test_dust()

    m.set_cylindrical_polar_grid_auto(399, 199, 1)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)


def test_manual_grid(tmpdir):
    """
    Regression test for the case where the analytical YSO grid is set manually.
    """

    dust = get_test_dust()

    m = AnalyticalYSOModel()

    m.star.radius = 1.
    m.star.temperature = 1000.
    m.star.luminosity = 1.

    d = m.add_flared_disk()
    d.mass = 1.
    d.rmin = 0.1
    d.rmax = 10.
    d.p = -1
    d.beta = 1.25
    d.h_0 = 0.1
    d.r_0 = 2.
    d.dust = dust

    r = np.linspace(0., 10., 100)
    t = [0., np.pi]
    p = [0., 2. * np.pi]

    m.set_spherical_polar_grid(r, t, p)

    m.set_n_photons(initial=0, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)