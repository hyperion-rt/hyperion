from astropy.tests.helper import pytest
from astropy import units as u

from numpy.testing import assert_allclose

from .. import ImageConf
from ...util.functions import virtual_file


class TestImageConf(object):

    def test_valid(self):
        i = ImageConf()
        i.set_image_size(512, 512)
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_wavelength_range(250, 1., 1000.)
        i.write(virtual_file())

    def test_missing_size(self):
        i = ImageConf()
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_wavelength_range(250, 1., 1000.)
        with pytest.raises(Exception) as exc:
            i.write(virtual_file())
        assert exc.value.args[0] == "Image size has not been set"

    def test_missing_limits(self):
        i = ImageConf()
        i.set_image_size(512, 512)
        i.set_wavelength_range(250, 1., 1000.)
        with pytest.raises(Exception) as exc:
            i.write(virtual_file())
        assert exc.value.args[0] == "Image limits have not been set"

    def test_missing_wavelength(self):
        i = ImageConf()
        i.set_image_size(512, 512)
        i.set_image_limits(-1., 1., -1., 1.)
        with pytest.raises(Exception) as exc:
            i.write(virtual_file())
        assert exc.value.args[0] == "Wavelength range has not been set"

    def test_filters(self):
        f = virtual_file()
        i = ImageConf()
        i.set_image_size(512, 512)
        i.set_image_limits(-1., 1., -1., 1.)

        f1 = i.add_filter()
        f1.name = 'F1'
        f1.spectral_coord = [1, 1.1, 1.2, 1.3] * u.micron
        f1.transmission = [0., 100., 50, 0.] * u.percent
        f1.detector_type = 'photons'
        f1.alpha = 0.
        f1.central_spectral_coord = 1.15 * u.micron

        f2 = i.add_filter()
        f2.name = 'F2'
        f2.spectral_coord = [2, 2.1, 2.2, 2.3, 2.4] * u.micron
        f2.transmission = [0., 50, 100, 60, 0.] * u.percent
        f2.detector_type = 'energy'
        f2.alpha = 1.
        f2.central_spectral_coord = 2.15 * u.micron

        i.write(f)
        i2 = ImageConf.read(f)

        # The [::-1] below are because the filters get sorted by frequency when
        # written out.

        assert_allclose(i2._filters[0].spectral_coord.to(u.Hz, equivalencies=u.spectral()).value,
                        i._filters[0].spectral_coord.to(u.Hz, equivalencies=u.spectral()).value[::-1])

        assert_allclose(i2._filters[0].transmission.to(u.one).value,
                        i._filters[0].transmission.to(u.one).value[::-1])

        assert_allclose(i2._filters[1].spectral_coord.to(u.Hz, equivalencies=u.spectral()).value,
                        i._filters[1].spectral_coord.to(u.Hz, equivalencies=u.spectral()).value[::-1])

        assert_allclose(i2._filters[1].transmission.to(u.one).value,
                        i._filters[1].transmission.to(u.one).value[::-1])

        with pytest.raises(ValueError) as exc:
            i.set_wavelength_range(250, 1., 1000.)
            i.write(f)
        assert exc.value.args[0] == "Cannot specify both filters and wavelength range"
