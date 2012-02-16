import pytest

from hyperion.densities import FlaredDisk, PowerLawEnvelope
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
    d = PowerLawEnvelope()
    d.__setattr__(parameter, 1.)


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_negative(parameter):
    d = PowerLawEnvelope()
    if parameter in ['power']:
        d.__setattr__(parameter, -1.)
    else:
        with pytest.raises(ValueError):
            d.__setattr__(parameter, -1.)  # negative values are not valid


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_optthin(parameter):
    d = PowerLawEnvelope()
    if parameter in ['rmin', 'rmax']:
        d.__setattr__(parameter, OptThinRadius(1.))
    else:
        with pytest.raises(ValueError):
            d.__setattr__(parameter, OptThinRadius(1.))  # not valid for these parameters


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid1(parameter):
    d = PowerLawEnvelope()
    with pytest.raises(ValueError):
        d.__setattr__(parameter, 'a')  # can't be string


@pytest.mark.parametrize(('parameter'), ['mass', 'rmin', 'rmax', 'power', 'rho_0', 'r_0'])
def test_power_law_envelope_invalid2(parameter):
    d = PowerLawEnvelope()
    with pytest.raises(ValueError):
        d.__setattr__(parameter, [1., 2.])  # should be scalar
