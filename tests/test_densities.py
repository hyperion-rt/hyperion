import pytest
import numpy as np

from hyperion.densities import FlaredDisk, PowerLawEnvelope, UlrichEnvelope, BipolarCavity
from hyperion.util.convenience import OptThinRadius

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
        with pytest.raises(ValueError):
            d.__setattr__(parameter, -1.)  # negative values are not valid


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_optthin(parameter):
    d = FlaredDisk()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError):
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_invalid1(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError):
        d.__setattr__(parameter, 'a')  # can't be string


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'p', 'beta', 'h_0', 'r_0'])
def test_flared_disk_invalid2(parameter):
    d = FlaredDisk()
    with pytest.raises(ValueError):
        d.__setattr__(parameter, [1., 2.])  # should be scalar

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
        with pytest.raises(ValueError):
            e.__setattr__(parameter, -1.)  # negative values are not valid


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_optthin(parameter):
    e = PowerLawEnvelope()
    if parameter in ['rmin', 'rmax']:
        e.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError):
            e.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid1(parameter):
    e = PowerLawEnvelope()
    with pytest.raises(ValueError):
        e.__setattr__(parameter, 'a')  # can't be string


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid2(parameter):
    e = PowerLawEnvelope()
    with pytest.raises(ValueError):
        e.__setattr__(parameter, [1., 2.])  # should be scalar


def test_power_law_envelope_swap1():
    e = PowerLawEnvelope()
    e.mass = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mass')
    e.rho_0 = 1.
    assert not hasattr(e, 'mass') and hasattr(e, 'rho_0')
    e.mass = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mass')


@pytest.mark.xfail
def test_power_law_envelope_swap2():
    e = PowerLawEnvelope()
    e.mass = 0.
    assert e.rho_0 == 0.


@pytest.mark.xfail
def test_power_law_envelope_swap3():
    e = PowerLawEnvelope()
    e.rho_0 = 0.
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
    with pytest.raises(ValueError):
        e.cavity = 1.  # should be BipolarCavity instance


def test_power_law_cavity_invalid2():
    e = PowerLawEnvelope()
    with pytest.raises(ValueError):
        e.cavity = np.array([1, 2, 3])  # should be BipolarCavity instance


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
        with pytest.raises(ValueError):
            e.__setattr__(parameter, -1.)  # negative values are not valid


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_optthin(parameter):
    e = UlrichEnvelope()
    if parameter in ['rmin', 'rmax']:
        e.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError):
            e.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_invalid1(parameter):
    e = UlrichEnvelope()
    with pytest.raises(ValueError):
        e.__setattr__(parameter, 'a')  # can't be string


@pytest.mark.parametrize(('parameter'), ['mdot', 'rmin', 'rmax', 'rho_0', 'rc', 'rho_amb'])
def test_ulrich_envelope_invalid2(parameter):
    e = UlrichEnvelope()
    with pytest.raises(ValueError):
        e.__setattr__(parameter, [1., 2.])  # should be scalar


def test_ulrich_envelope_swap():
    e = UlrichEnvelope()
    e.mdot = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mdot')
    e.rho_0 = 1.
    assert not hasattr(e, 'mdot') and hasattr(e, 'rho_0')
    e.mdot = 1.
    assert not hasattr(e, 'rho_0') and hasattr(e, 'mdot')


@pytest.mark.xfail
def test_ulrich_envelope_swap2():
    e = UlrichEnvelope()
    e.mdot = 0.
    assert e.rho_0 == 0.


@pytest.mark.xfail
def test_ulrich_envelope_swap3():
    e = UlrichEnvelope()
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
    with pytest.raises(ValueError):
        e.cavity = 1.  # should be BipolarCavity instance


def test_ulrich_cavity_invalid2():
    e = UlrichEnvelope()
    with pytest.raises(ValueError):
        e.cavity = np.array([1, 2, 3])  # should be BipolarCavity instance
