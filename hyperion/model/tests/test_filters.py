
import os
import shutil
import tempfile

import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from astropy import units as u

from .. import Model
from ..image import Image
from ...util.functions import random_id
from .test_helpers import get_test_dust, get_highly_reflective_dust


class TestFilters(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        m.add_density_grid(np.array([[[1.]]]), get_test_dust())

        s = m.add_point_source()
        s.name = 'first'
        s.luminosity = 1.
        s.temperature = 6000.

        s = m.add_point_source()
        s.name = 'second'
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=True, image=True)
        i.set_viewing_angles([1., 2., 3.], [1., 2., 3.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)

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

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1000)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    @pytest.mark.requires_hyperion_binaries
    def test_image_wav(self):
        image = self.m.get_image()
        np.testing.assert_allclose(image.nu, [2.60689094e+14, 1.39438353e+14])
        np.testing.assert_allclose(image.wav, [1.15, 2.15])

    @pytest.mark.requires_hyperion_binaries
    def test_sed_wav(self):
        sed = self.m.get_sed()
        np.testing.assert_allclose(sed.nu, [2.60689094e+14, 1.39438353e+14])
        np.testing.assert_allclose(sed.wav, [1.15, 2.15])

    @pytest.mark.requires_hyperion_binaries
    def test_image_shape(self):
        image = self.m.get_image()
        assert image.val.shape == (3, 20, 10, 2)

    @pytest.mark.requires_hyperion_binaries
    def test_sed_shape(self):
        sed = self.m.get_sed()
        assert sed.val.shape == (3, 1, 2)

    @pytest.mark.requires_hyperion_binaries
    def test_image_values(self):
        image = self.m.get_image(units='MJy/sr', distance=1)
        np.testing.assert_allclose(np.sum(image.val[:, :, :, 0]), 3438.059082285024, rtol=0.1)
        np.testing.assert_allclose(np.sum(image.val[:, :, :, 1]), 2396.4803378036186, rtol=0.1)

    @pytest.mark.requires_hyperion_binaries
    def test_image_filter_names(self):
        image = self.m.get_image()
        assert image.filter_names == ['F1', 'F2']

    @pytest.mark.requires_hyperion_binaries
    def test_sed_filter_names(self):
        sed = self.m.get_sed()
        assert sed.filter_names == ['F1', 'F2']


@pytest.mark.requires_hyperion_binaries
def test_filter_names_linked_input(tmpdir):

    # If the input is linked rather than copied into the output file, the
    # filter names should be read through the link, and should be None if the
    # input file is no longer available

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), get_test_dust())

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([1.], [1.])

    f = i.add_filter()
    f.name = 'F1'
    f.spectral_coord = [1, 1.1, 1.2, 1.3] * u.micron
    f.transmission = [0., 100., 50, 0.] * u.percent
    f.detector_type = 'photons'
    f.alpha = 0.
    f.central_spectral_coord = 1.15 * u.micron

    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=1000)
    m.set_copy_input(False)

    input_file = tmpdir.join('test.rtin').strpath
    m.write(input_file)
    mo = m.run(tmpdir.join('test.rtout').strpath)

    assert mo.get_sed().filter_names == ['F1']

    os.rename(input_file, tmpdir.join('moved.rtin').strpath)

    sed = mo.get_sed()
    assert sed.filter_names is None
    np.testing.assert_allclose(sed.wav, [1.15])
