from __future__ import print_function, division

import pytest
import numpy as np

from .. import FlaredDisk, PowerLawEnvelope, UlrichEnvelope, BipolarCavity
from ...util.convenience import OptThinRadius

# A fake star class so that star.mass is defined
class Star(object):
    def __init__(self):
        self.mass = None

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
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mass')
    e.rho_0 = 1.
    assert not hasattr(e, 'mass') and hasattr(e, 'rho_0')
    e.mass = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mass')


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


def test_power_law_cavity():
    e = PowerLawEnvelope()
    e.add_bipolar_cavity()
    assert e.cavity.envelope is e


def test_power_law_cavity_direct():
    e = PowerLawEnvelope()
    e.cavity = BipolarCavity()
    assert e.cavity.envelope is e


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


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_positive(parameter):
    e = UlrichEnvelope()
    e.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_negative(parameter):
    e = UlrichEnvelope()
    if parameter in ['power']:
        e.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, -1.)  # negative values are not valid
        assert exc.value.args[0] == parameter + ' should be positive'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_optthin(parameter):
    e = UlrichEnvelope()
    if parameter in ['rmin', 'rmax']:
        e.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError) as exc:
            e.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters
        assert exc.value.args[0] == parameter + ' should be a scalar value'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_invalid1(parameter):
    e = UlrichEnvelope()
    with pytest.raises(ValueError) as exc:
        e.__setattr__(parameter, 'a')  # can't be string
    if parameter in ['rmin', 'rmax']:
        assert exc.value.args[0] == parameter + ' should be a numerical value or an OptThinRadius instance'
    else:
        assert exc.value.args[0] == parameter + ' should be a numerical value'


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
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
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mdot')
    e.rho_0 = 1.
    assert not hasattr(e, 'mdot') and hasattr(e, 'rho_0')
    e.mdot = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mdot')


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


def test_ulrich_cavity():
    e = UlrichEnvelope()
    e.add_bipolar_cavity()
    assert e.cavity.envelope is e


def test_ulrich_cavity_direct():
    e = UlrichEnvelope()
    e.cavity = BipolarCavity()
    assert e.cavity.envelope is e


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
