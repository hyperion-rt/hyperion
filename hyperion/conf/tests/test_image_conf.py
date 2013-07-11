from astropy.tests.helper import pytest

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
