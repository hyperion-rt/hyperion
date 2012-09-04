from __future__ import print_function, division

import pytest
import numpy as np

from .. import AnalyticalYSOModel
from ...util.constants import msun, rsun, lsun, tsun, au
from ...util.functions import random_filename
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


def test_analytical_yso_full():

    m = basic_analytical_model()

    m.set_spherical_polar_grid_auto(10, 10, 10)

    m.set_n_photons(initial=100, imaging=100)

    m.write(random_filename())


def test_analytical_yso_nogrid_invalid():

    m = AnalyticalYSOModel()
    with pytest.raises(Exception) as e:
        m.write(random_filename())
    assert e.value.args[0] == 'The coordinate grid needs to be defined before calling AnalyticalModelYSO.write(...)'


def test_analytical_yso_nostar_invalid():

    m = AnalyticalYSOModel()
    with pytest.raises(Exception) as e:
        m.set_spherical_polar_grid_auto(1, 1, 1)
    assert e.value.args[0] == 'The central source radius need to be defined before the grid can be set up'


def test_analytical_yso_optthinradius():

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

    m.write(random_filename())


def test_analytical_yso_optthinradius_check_frozen():

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

    with pytest.raises(Exception) as exc:
        e.mass = 1.
    assert "Attribute mass can no longer be changed"

    with pytest.raises(Exception) as exc:
        e.dust = 'test'
    assert "Attribute dust can no longer be changed"

    with pytest.raises(Exception) as exc:
        m.star.radius = 1.
    assert "Attribute radius can no longer be changed"

    m.set_n_photons(initial=100, imaging=100)

    m.write(random_filename())


def test_analytical_yso_optthinradius_manual_check_frozen():

    dust = get_test_dust()

    m = basic_analytical_model()

    e = m.add_power_law_envelope()
    e.mass = 1.e-2 * msun
    e.rmin = OptThinRadius(1000.)
    e.r_0 = au
    e.rmax = OptThinRadius(10.)
    e.power = -2.
    e.dust = dust

    assert not np.isscalar(e.rmin)
    assert not np.isscalar(e.rmax)

    m.evaluate_optically_thin_radii()

    assert np.isscalar(e.rmin)
    assert np.isscalar(e.rmax)

    with pytest.raises(Exception) as exc:
        e.mass = 1.
    assert "Attribute mass can no longer be changed"

    with pytest.raises(Exception) as exc:
        e.dust = 'test'
    assert "Attribute dust can no longer be changed"

    with pytest.raises(Exception) as exc:
        m.star.radius = 1.
    assert "Attribute radius can no longer be changed"

    m.set_spherical_polar_grid([0., 0.1 * au, 10. * au], [0., np.pi], [0., 2. * np.pi])

    m.set_n_photons(initial=100, imaging=100)

    m.write(random_filename())


def test_analytical_yso_add_density():
    m = AnalyticalYSOModel()
    with pytest.raises(NotImplementedError) as exc:
        m.add_density_grid()
    assert exc.value.args[0] == 'add_density_grid cannot be used for AnalyticalYSOModel'


def test_analytical_yso_use_quantities_invalid():

    output_file = random_filename()

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())
    m.run(output_file)

    m2 = basic_analytical_model()
    m2.set_spherical_polar_grid_auto(10, 10, 10)
    m2.set_n_photons(initial=100, imaging=100)
    with pytest.raises(NotImplementedError) as exc:
        m2.use_quantities(output_file)
    assert exc.value.args[0] == "Cannot use previous density in AnalyticalYSOModel. If you want to use just the previous specific_energy, specify quantities=['specific_energy']."


def test_analytical_yso_use_quantities():

    output_file = random_filename()

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())
    m.run(output_file)

    m2 = basic_analytical_model()
    m2.set_spherical_polar_grid_auto(10, 10, 10)
    m2.set_n_photons(initial=100, imaging=100)
    m2.use_quantities(output_file, quantities=['specific_energy'])
    m2.write(random_filename())
    m2.run(random_filename())


def test_analytical_yso_use_geometry():

    output_file = random_filename()

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())
    m.run(output_file)

    m2 = basic_analytical_model()
    m2.use_geometry(output_file)
    m2.set_n_photons(initial=100, imaging=100)
    m2.write(random_filename())
    m2.run(random_filename())


def test_analytical_yso_use_geometry_quantities():

    output_file = random_filename()

    m = basic_analytical_model()
    m.set_spherical_polar_grid_auto(10, 10, 10)
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())
    m.run(output_file)

    m2 = basic_analytical_model()
    m2.use_geometry(output_file)
    m2.set_n_photons(initial=100, imaging=100)
    m2.use_quantities(output_file, quantities=['specific_energy'])
    m2.write(random_filename())
    m2.run(random_filename())


def test_ambient_medium():

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

    m.write(random_filename())


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
        with pytest.raises(ValueError) as exc:
            m.set_cylindrical_polar_grid_auto(100, 20, 3)
        assert exc.value.args[0] == "R_min is 0, so cannot set up the grid cell walls automatically. Use set_cylindrical_polar_grid() instead to specify the cell wall positions."
    else:
        with pytest.raises(ValueError) as exc:
            m.set_spherical_polar_grid_auto(100, 20, 3)
        assert exc.value.args[0] == "R_min is 0, so cannot set up the grid cell walls automatically. Use set_spherical_polar_grid() instead to specify the cell wall positions."


def test_complete():

    m = AnalyticalYSOModel()

    m.star.radius = 1.
    m.star.temperature = 1000.
    m.star.luminosity = 1.

    d = m.add_flared_disk()
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

    m.write(random_filename())
