from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from .. import FlaredDisk, AlphaDisk, PowerLawEnvelope, UlrichEnvelope, BipolarCavity, AmbientMedium
from ...grid import SphericalPolarGrid
from ...util.convenience import OptThinRadius
from ...util.constants import G

# A fake star class so that star.mass is defined


class Star(object):
    def __init__(self):
        self.mass = None
        self.radius = None

# Flared Disk


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0'])
def test_flared_disk_positive(parameter):
    d = FlaredDisk()
    d.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0'])
def test_flared_disk_negative(parameter):
    d = FlaredDisk()
    if parameter in ['p', 'beta']:
        d.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0'])
def test_flared_disk_optthin(parameter):
    d = FlaredDisk()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0'])
def test_flared_disk_invalid1(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0'])
def test_flared_disk_invalid2(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


def test_flared_disk_mass_swap1():
    e = FlaredDisk()
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None
    e.rho_0 = 1.
    assert e._rho_0 is not None and e._mass is None
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None


def test_flared_disk_mass_swap2():
    e = FlaredDisk()
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mass = 0.
    assert e.rho_0 == 0.


def test_flared_disk_mass_swap3():
    e = FlaredDisk()
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.rho_0 = 0.
    assert e.mass == 0.


def test_flared_disk_mass_swap_invertible():
    e = FlaredDisk()
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mass = 4.792849
    rho_0 = e.rho_0
    e.mass = 0.
    e.rho_0 = rho_0
    assert e.mass == 4.792849

# Alpha Disk


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_positive(parameter):
    d = AlphaDisk()
    d.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_negative(parameter):
    d = AlphaDisk()
    if parameter in ['p', 'beta']:
        d.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_optthin(parameter):
    d = AlphaDisk()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_invalid1(parameter):
    d = AlphaDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'rho_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_invalid2(parameter):
    d = AlphaDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


def test_alpha_disk_mass_swap1():
    e = AlphaDisk()
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None
    e.rho_0 = 1.
    assert e._rho_0 is not None and e._mass is None
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None


def test_alpha_disk_mass_swap2():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 1.
    e.star.radius = 1.
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mass = 0.
    e.lvisc = 0.
    assert e.rho_0 == 0.


def test_alpha_disk_mass_swap3():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 1.
    e.star.radius = 1.
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.rho_0 = 0.
    e.lvisc = 0.
    assert e.mass == 0.


def test_alpha_disk_mass_swap_invertible():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 1.
    e.star.radius = 1.
    e.rmin = 1.
    e.rmax = 10.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.lvisc = 0.
    e.mass = 4.792849
    rho_0 = e.rho_0
    e.mass = 0.
    e.rho_0 = rho_0
    assert e.mass == 4.792849


def test_alpha_disk_swap1():
    e = AlphaDisk()
    e.mdot = 1.
    assert e._lvisc is None and e._mdot is not None
    e.lvisc = 1.
    assert e._lvisc is not None and e._mdot is None
    e.mdot = 1.
    assert e._lvisc is None and e._mdot is not None


def test_alpha_disk_swap2():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 1.
    e.star.radius = 1.
    e.rmin = 1.
    e.rmax = 10.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mdot = 0.
    assert e.lvisc == 0.


def test_alpha_disk_swap3():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 1.
    e.star.radius = 1.
    e.rmin = 1.
    e.rmax = 10.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.lvisc = 0.
    assert e.mdot == 0.


def test_alpha_disk_lvisc_calc():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 0.5
    e.star.radius = 0.5
    e.rmin = 1.
    e.rmax = 2.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mdot = 4.
    assert_array_almost_equal_nulp(e.lvisc, G * (1.5 - 2. / np.sqrt(2.) + 0.5), 2)


def test_alpha_disk_mdot_calc():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 0.5
    e.star.radius = 0.5
    e.rmin = 1.
    e.rmax = 2.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.lvisc = G * (1.5 - 2. / np.sqrt(2.) + 0.5)
    assert_array_almost_equal_nulp(e.mdot, 4., 2)


def test_alpha_disk_lvisc_map():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 0.5
    e.star.radius = 0.5
    e.rmin = 1.
    e.rmax = 2.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.lvisc = G * (1.5 - 2. / np.sqrt(2.) + 0.5)

    # Set up grid
    r = np.linspace(0., e.rmax * 2., 200)
    t = np.linspace(0., np.pi, 200)
    p = np.linspace(0., 2. * np.pi, 5)
    g = SphericalPolarGrid(r, t, p)

    ref = G * (1.5 - 2. / np.sqrt(2.) + 0.5)
    act = np.sum(e.accretion_luminosity(g))

    assert abs(act - ref) / ref < 1.e-3


def test_alpha_disk_mdot_map():
    e = AlphaDisk()
    e.star = Star()
    e.star.mass = 0.5
    e.star.radius = 0.5
    e.rmin = 1.
    e.rmax = 2.
    e.mass = 1.
    e.r_0 = 5.
    e.h_0 = 1.
    e.p = -1.
    e.beta = 1.25
    e.mdot = 4.

    # Set up grid
    r = np.linspace(0., e.rmax * 2., 200)
    t = np.linspace(0., np.pi, 200)
    p = np.linspace(0., 2. * np.pi, 5)
    g = SphericalPolarGrid(r, t, p)

    ref = G * (1.5 - 2. / np.sqrt(2.) + 0.5)
    act = np.sum(e.accretion_luminosity(g))

    assert abs(act - ref) / ref < 1.e-3

# Power Law Envelope


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_positive(parameter):
    e = PowerLawEnvelope()
    e.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_negative(parameter):
    e = PowerLawEnvelope()
    if parameter in ['power']:
        e.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_optthin(parameter):
    e = PowerLawEnvelope()
    if parameter in ['rmin', 'rmax']:
        e.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid1(parameter):
    e = PowerLawEnvelope()
    with pytest.raises(ValueError) as exc:
        e.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid2(parameter):
    e = PowerLawEnvelope()
    with pytest.raises(ValueError) as exc:
        e.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


def test_power_law_envelope_swap1():
    e = PowerLawEnvelope()
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None
    e.rho_0 = 1.
    assert e._rho_0 is not None and e._mass is None
    e.mass = 1.
    assert e._rho_0 is None and e._mass is not None


def test_power_law_envelope_swap2():
    e = PowerLawEnvelope()
    e.mass = 0.
    e.rmin = 1.
    e.r_0 = 2.
    e.rmax = 10.
    e.power = -2.
    assert e.rho_0 == 0.


def test_power_law_envelope_swap3():
    e = PowerLawEnvelope()
    e.rho_0 = 0.
    e.rmin = 1.
    e.r_0 = 2.
    e.rmax = 10.
    e.power = -2.
    assert e.mass == 0.


def test_power_law_envelope_rho0_calc():
    e = PowerLawEnvelope()
    e.rmin = 9.
    e.rmax = 18.
    e.r_0 = 3.
    e.power = -1.
    e.mass = 4. * np.pi * 27 * 27
    assert_array_almost_equal_nulp(e.rho_0, 2., 2)


def test_power_law_envelope_mass_calc():
    e = PowerLawEnvelope()
    e.rmin = 9.
    e.rmax = 18.
    e.r_0 = 3.
    e.power = -1.
    e.rho_0 = 2.
    assert_array_almost_equal_nulp(e.mass, 4. * np.pi * 27 * 27, 2)


def test_power_law_cavity():
    e = PowerLawEnvelope()
    e.add_bipolar_cavity()
    assert e.cavity._envelope is e


def test_power_law_cavity_direct():
    e = PowerLawEnvelope()
    e.cavity = BipolarCavity()
    assert e.cavity._envelope is e


def test_power_law_cavity_invalid1():
    e = PowerLawEnvelope()
    with pytest.raises(ValueError) as exc:
        e.cavity = 1.  # should be BipolarCavity instance
    assert exc.value.args[0] == 'cavity should be an instance of BipolarCavity'


def test_power_law_cavity_invalid2():
    e = PowerLawEnvelope()
    with pytest.raises(ValueError) as exc:
        e.cavity = np.array([1, 2, 3])  # should be BipolarCavity instance
    assert exc.value.args[0] == 'cavity should be an instance of BipolarCavity'


# Ulrich Envelope


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc'])
def test_ulrich_envelope_positive(parameter):
    e = UlrichEnvelope()
    e.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc'])
def test_ulrich_envelope_negative(parameter):
    e = UlrichEnvelope()
    if parameter in ['power']:
        e.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc'])
def test_ulrich_envelope_optthin(parameter):
    e = UlrichEnvelope()
    if parameter in ['rmin', 'rmax']:
        e.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc'])
def test_ulrich_envelope_invalid1(parameter):
    e = UlrichEnvelope()
    with pytest.raises(ValueError) as exc:
        e.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc'])
def test_ulrich_envelope_invalid2(parameter):
    e = UlrichEnvelope()
    with pytest.raises(ValueError) as exc:
        e.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


def test_ulrich_envelope_swap():
    e = UlrichEnvelope()
    e.mdot = 1.
    assert e._rho_0 is None and e._mdot is not None
    e.rho_0 = 1.
    assert e._rho_0 is not None and e._mdot is None
    e.mdot = 1.
    assert e._rho_0 is None and e._mdot is not None


def test_ulrich_envelope_swap2():
    e = UlrichEnvelope()
    e.star = Star()
    e.star.mass = 1.
    e.rmin = 1.
    e.rc = 2.
    e.rmax = 10.
    e.mdot = 0.
    assert e.rho_0 == 0.


def test_ulrich_envelope_swap3():
    e = UlrichEnvelope()
    e.star = Star()
    e.star.mass = 1.
    e.rmin = 1.
    e.rc = 2.
    e.rmax = 10.
    e.rho_0 = 0.
    assert e.mdot == 0.


def test_ulrich_envelope_rho0_calc():
    e = UlrichEnvelope()
    e.star = Star()
    e.star.mass = 10.
    e.rmin = 1.
    e.rc = 2.
    e.rmax = 10.
    e.mdot = 7.
    assert_array_almost_equal_nulp(e.rho_0, 7. / (4. * np.pi * np.sqrt(80. * G)), 2)


def test_ulrich_envelope_mdot_calc():
    e = UlrichEnvelope()
    e.star = Star()
    e.star.mass = 10.
    e.rmin = 1.
    e.rc = 2.
    e.rmax = 10.
    e.rho_0 = 7. / (4. * np.pi * np.sqrt(80. * G))
    assert_array_almost_equal_nulp(e.mdot, 7., 2)


def test_ulrich_cavity():
    e = UlrichEnvelope()
    e.add_bipolar_cavity()
    assert e.cavity._envelope is e


def test_ulrich_cavity_direct():
    e = UlrichEnvelope()
    e.cavity = BipolarCavity()
    assert e.cavity._envelope is e


def test_ulrich_cavity_invalid1():
    e = UlrichEnvelope()
    with pytest.raises(ValueError) as exc:
        e.cavity = 1.  # should be BipolarCavity instance
    assert exc.value.args[0] == 'cavity should be an instance of BipolarCavity'


def test_ulrich_cavity_invalid2():
    e = UlrichEnvelope()
    with pytest.raises(ValueError) as exc:
        e.cavity = np.array([1, 2, 3])  # should be BipolarCavity instance
    assert exc.value.args[0] == 'cavity should be an instance of BipolarCavity'


# Bipolar Cavities

@pytest.mark.parametrize(('parameter'), ['theta_0', 'r_0', 'rho_0', 'rho_exp'])
def test_bipolar_cavity_positive(parameter):
    c = BipolarCavity()
    c.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['theta_0', 'r_0', 'rho_0', 'rho_exp'])
def test_bipolar_cavity_negative(parameter):
    c = BipolarCavity()
    if parameter in ['power', 'rho_exp']:
        c.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            c.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['theta_0', 'r_0', 'rho_0', 'rho_exp'])
def test_bipolar_cavity_invalid1(parameter):
    c = BipolarCavity()
    with pytest.raises(ValueError) as exc:
        c.__setattr__(parameter, 'a')  # can't be string
    assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['theta_0', 'r_0', 'rho_0', 'rho_exp'])
def test_bipolar_cavity_invalid2(parameter):
    c = BipolarCavity()
    with pytest.raises(ValueError) as exc:
        c.__setattr__(parameter, [1., 2.])  # should be scalar
    assert exc.value.args[0] == parameter + ' should be a scalar value'


# Ambient medium


@pytest.mark.parametrize(('parameter'), ['rho', 'rmin', 'rmax'])
def test_ambient_positive(parameter):
    a = AmbientMedium()
    a.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['rho', 'rmin', 'rmax'])
def test_ambient_negative(parameter):
    a = AmbientMedium()
    with pytest.raises(ValueError) as exc:
        a.__setattr__(parameter, -1.)  # negative values are not valid
    assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['rho', 'rmin', 'rmax'])
def test_ambient_optthin(parameter):
    a = AmbientMedium()
    if parameter in ['rmin', 'rmax']:
        a.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            a.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['rho', 'rmin', 'rmax'])
def test_ambient_invalid1(parameter):
    a = AmbientMedium()
    with pytest.raises(ValueError) as exc:
        a.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['rho', 'rmin', 'rmax'])
def test_ambient_invalid2(parameter):
    a = AmbientMedium()
    with pytest.raises(ValueError) as exc:
        a.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


def test_ambient_densities_1():

    r = np.linspace(0., 10., 10)
    t = [0., np.pi]
    p = [0., 2 * np.pi]
    g = SphericalPolarGrid(r, t, p)

    a = AmbientMedium()
    a.rho = 2.
    a.rmin = 0.1
    a.rmax = 10.

    expected = np.ones(len(r) - 1) * 2.
    assert_array_almost_equal_nulp(a.density(g)[0,0,:], expected, 10)

def test_ambient_densities_2():

    r = np.linspace(0., 10., 10)
    t = [0., np.pi]
    p = [0., 2 * np.pi]
    g = SphericalPolarGrid(r, t, p)

    # Set up envelope
    p1 = PowerLawEnvelope()
    p1.power = -2
    p1.r_0 = 1.
    p1.rho_0 = 10.
    p1.rmin = 0.1
    p1.rmax = 10.

    a = AmbientMedium()
    a.rho = 2.
    a.rmin = 0.1
    a.rmax = 10.
    a.subtract = [p1]

    expected = 10. * g.r ** -2
    assert_array_almost_equal_nulp(p1.density(g)[0,0,:], expected, 10)

    expected = np.clip(2 - 10. * g.r ** -2, 0., np.inf)
    assert_array_almost_equal_nulp(a.density(g)[0,0,:], expected, 10)

    expected = np.clip(10. * g.r ** -2, 2., np.inf)
    assert_array_almost_equal_nulp((a.density(g) + p1.density(g))[0,0,:],
                                   expected, 10)

def test_ambient_densities_3():

    r = np.linspace(0., 10., 10)
    t = [0., np.pi]
    p = [0., 2 * np.pi]
    g = SphericalPolarGrid(r, t, p)

    # Set up envelope
    p1 = PowerLawEnvelope()
    p1.power = -2
    p1.r_0 = 1.
    p1.rho_0 = 10.
    p1.rmin = 0.1
    p1.rmax = 10.

    # Set up another envelope
    p2 = PowerLawEnvelope()
    p2.power = -1.5
    p2.r_0 = 1.
    p2.rho_0 = 8.
    p2.rmin = 0.1
    p2.rmax = 10.

    a = AmbientMedium()
    a.rho = 2.
    a.rmin = 0.1
    a.rmax = 10.
    a.subtract = [p1, p2]

    expected = 10. * g.r ** -2
    assert_array_almost_equal_nulp(p1.density(g)[0,0,:], expected, 10)

    expected = np.clip(2 - 10. * g.r ** -2 - 8. * g.r ** -1.5, 0., np.inf )
    assert_array_almost_equal_nulp(a.density(g)[0,0,:], expected, 10)

    expected = np.clip(10. * g.r ** -2 + 8. * g.r ** -1.5, 2., np.inf)
    assert_array_almost_equal_nulp((a.density(g) + p1.density(g) + p2.density(g))[0,0,:],
                                   expected, 10)


def test_ambient_densities_4():

    # Regression test for #106 - could not add a bipolar cavity to the list of
    # densities to subtract from the ambient medium.

    r = np.linspace(0., 10., 10)
    t = np.linspace(0., np.pi,10)
    p = [0., 2 * np.pi]
    g = SphericalPolarGrid(r, t, p)

    # Set up envelope
    p1 = PowerLawEnvelope()
    p1.power = -2
    p1.r_0 = 1.
    p1.rho_0 = 10.
    p1.rmin = 0.1
    p1.rmax = 10.

    # Set up bipolar cavity
    c = p1.add_bipolar_cavity()
    c.power = 2
    c.rho_0 = 3.
    c.r_0 = 1.
    c.theta_0 = 80.

    a = AmbientMedium()
    a.rho = 2.
    a.rmin = 0.1
    a.rmax = 10.
    a.subtract = [p1, c]

    expected = np.repeat(3., 9)
    assert_array_almost_equal_nulp((p1.density(g) + c.density(g) + a.density(g))[0,0,:], expected, 10)

    c.rho_0 = 2.

    expected = np.repeat(2., 9)
    assert_array_almost_equal_nulp((p1.density(g) + c.density(g) + a.density(g))[0,0,:], expected, 10)

    c.rho_0 = 1.

    expected = np.repeat(2., 9)
    assert_array_almost_equal_nulp((p1.density(g) + c.density(g) + a.density(g))[0,0,:], expected, 10)
