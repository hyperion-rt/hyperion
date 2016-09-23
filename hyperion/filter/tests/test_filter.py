import os

import pytest
import numpy as np

from astropy import units as u

from ...util.functions import virtual_file

from .. import Filter

ROOT = os.path.dirname(__file__)


def test_init():
    Filter()


def test_basic():

    f = Filter()

    f.name = "2J"
    f.spectral_coord = [1, 2, 3] * u.Hz
    f.transmission = [4, 5, 6] * u.percent
    f.central_spectral_coord = 1.5 * u.Hz
    f.detector_type = 'photons'
    f.alpha = -1.

    assert f.name == "2J"
    np.testing.assert_allclose(f.spectral_coord.to(u.Hz).value, [1., 2., 3.])
    np.testing.assert_allclose(f.transmission.to(u.percent).value, [4., 5., 6.])
    np.testing.assert_allclose(f.central_spectral_coord.to(u.Hz).value, 1.5)
    assert f.detector_type == 'photons'
    assert f.alpha == -1.


@pytest.mark.parametrize('value', [object(), 1., [1, 2, 3]])
def test_name_invalid_type(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.name = value
    assert exc.value.args[0] == 'name should be given as a string'


@pytest.mark.parametrize('value', [1. * u.Hz, [[1, 2, 3], [4, 5, 6]] * u.Hz])
def test_spectral_coord_invalid_shape(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.spectral_coord = value
    assert exc.value.args[0] == 'spectral_coord should be a 1-d sequence'


@pytest.mark.parametrize('value', ['string', object(), 1., [1, 2, 3]])
def test_spectral_coord_invalid_type(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.spectral_coord = value
    assert exc.value.args[0] == 'spectral_coord should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_spectral_coord_invalid_unit(unit):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.spectral_coord = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'spectral_coord should be given in units of frequency, length, energy'


@pytest.mark.parametrize('value', ['string', object(), 1., [1, 2, 3]])
def test_central_spectral_coord_invalid_type(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.central_spectral_coord = value
    assert exc.value.args[0] == 'spectral_coord should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_central_spectral_coord_invalid_unit(unit):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.central_spectral_coord = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'spectral_coord should be given in units of frequency, length, energy'


def test_roundtrip():

    f = Filter()

    f.name = "2J"
    f.spectral_coord = [1, 2, 3] * u.Hz
    f.transmission = [4, 5, 6] * u.percent
    f.central_spectral_coord = 1.5 * u.Hz
    f.detector_type = 'photons'
    f.alpha = -1.

    v = virtual_file()

    f.to_hdf5_group(v, "testing")
    f2 = Filter.from_hdf5_group(v, "testing")

    assert f.name == f2.name
    np.testing.assert_allclose(f.spectral_coord.to(u.Hz).value,
                               f2.spectral_coord.to(u.Hz).value)
    np.testing.assert_allclose(f.transmission.to(u.percent).value,
                               f2.transmission.to(u.percent).value)
