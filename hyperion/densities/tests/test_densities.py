from __future__ import print_function, division

import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from .. import FlaredDisk, AlphaDisk, PowerLawEnvelope, UlrichEnvelope, BipolarCavity
from ...grid import SphericalPolarGrid
from ...util.convenience import OptThinRadius
from ...util.constants import G

# A fake star class so that star.mass is defined


class Star(object):
    def __init__(self):
        self.mass = None
        self.radius = None

# Flared Disk


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_positive(parameter):
    d = FlaredDisk()
    d.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_negative(parameter):
    d = FlaredDisk()
    if parameter in ['p', 'beta']:
        d.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_optthin(parameter):
    d = FlaredDisk()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_invalid1(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_invalid2(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'

# Alpha Disk


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_positive(parameter):
    d = AlphaDisk()
    d.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_negative(parameter):
    d = AlphaDisk()
    if parameter in ['p', 'beta']:
        d.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_optthin(parameter):
    d = AlphaDisk()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_invalid1(parameter):
    d = AlphaDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0', 'mdot', 'lvisc'])
def test_alpha_disk_invalid2(parameter):
    d = AlphaDisk()
    with pytest.raises(ValueError) as exc:
        d.__setattr__(parameter, [1., 2.])  # should be scalar
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a scalar value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a scalar value'


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
