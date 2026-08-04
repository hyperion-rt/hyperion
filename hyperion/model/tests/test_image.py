from __future__ import print_function, division

import os
import shutil
import tempfile

import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from .. import Model
from ..image import Image
from ...util.functions import random_id
from .test_helpers import get_test_dust, get_highly_reflective_dust


@pytest.mark.requires_hyperion_binaries
class TestImageSimpleModel(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_stokes(True)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_image_group(self):
        wav, nufnu = self.m.get_image(group=0)

    def test_image_group_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(group=-2)
            # negative indexing allowed, but only one group present
        assert exc.value.args[0] == 'File only contains 1 image/SED group(s)'

    def test_image_group_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(group=1)
            # zero-based, and only one group present
        assert exc.value.args[0] == 'File only contains 1 image/SED group(s)'

    def test_image_dim(self):
        wav, nufnu = self.m.get_image()
        assert nufnu.shape == (2, 20, 10, 5)

    def test_image_dim_incl1(self):
        wav, nufnu = self.m.get_image(inclination=0)
        assert nufnu.shape == (20, 10, 5)

    def test_image_dim_incl2(self):
        wav, nufnu = self.m.get_image(inclination=1)
        assert nufnu.shape == (20, 10, 5)

    def test_image_dim_incl_invalid1(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_image(inclination=2)

    def test_image_dim_incl_invalid2(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_image(inclination=-3)

    def test_image_dim_incl_invalid3(self):
        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_image(inclination=12.3)
        assert exc.value.args[0] == "inclination should be an integer (it should be the index of the inclination, not the value itself)"

    @pytest.mark.parametrize(('stokes'), ['I', 'Q', 'U', 'V',
                                          'linpol', 'circpol'])
    def test_image_stokes(self, stokes):
        wav, nufnu = self.m.get_image(stokes=stokes)
        assert nufnu.shape == (2, 20, 10, 5)

    @pytest.mark.parametrize(('stokes'), ['A', 'b', 1, (3,),  # invalid values
                                          'i', 'q', 'u', 'v'])  # lowercase
    def test_image_stokes_invalid(self, stokes):
        with pytest.raises(ValueError):
            wav, nufnu = self.m.get_image(stokes=stokes)

    @pytest.mark.parametrize(('units'), ['ergs/s'])
    def test_image_nodistance_units(self, units):
        wav, nufnu = self.m.get_image(units=units)

    @pytest.mark.parametrize(('units'), ['ergs/cm^2/s', 'mJy', 'Jy', 'ergs/cm^2/s/Hz', 'MJy/sr'])
    def test_image_nodistance_units_invalid(self, units):
        with pytest.raises(ValueError):
            wav, nufnu = self.m.get_image(units=units)


@pytest.mark.requires_hyperion_binaries
class TestImageSimpleModelTrackingDetailed(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        m.add_density_grid(np.array([[[1.e-30]]]), get_test_dust())

        s = m.add_point_source()
        s.name = 'first'
        s.luminosity = 1.
        s.temperature = 6000.

        s = m.add_point_source()
        s.name = 'second'
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_track_origin('detailed')

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_image_source_all(self):
        wav, nufnu = self.m.get_image(source_id='all', component='source_emit')

    def test_image_source_valid1(self):
        wav, nufnu = self.m.get_image(source_id=0, component='source_emit')

    def test_image_source_valid2(self):
        wav, nufnu = self.m.get_image(source_id=1, component='source_emit')

    def test_image_source_valid3(self):
        wav, nufnu = self.m.get_image(source_id='first', component='source_emit')

    def test_image_source_valid4(self):
        wav, nufnu = self.m.get_image(source_id='second', component='source_emit')

    def test_image_source_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(source_id=-1, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

    def test_image_source_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(source_id=2, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

    def test_image_source_invalid3(self):
        with pytest.raises(ValueError) as exc:
            self.m.get_image(component='source')
        assert exc.value.args[0] == "component should be one of total/source_emit/dust_emit/source_scat/dust_scat since track_origin='detailed'"

    def test_image_dust_all(self):
        wav, nufnu = self.m.get_image(dust_id='all', component='dust_emit')

    def test_image_dust_valid1(self):
        wav, nufnu = self.m.get_image(dust_id=0, component='dust_emit')

    def test_image_dust_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(dust_id=-1, component='dust_emit')
        assert exc.value.args[0] == 'dust_id should be between 0 and 0'

    def test_image_dust_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(dust_id=1, component='dust_emit')
        assert exc.value.args[0] == 'dust_id should be between 0 and 0'


@pytest.mark.requires_hyperion_binaries
class TestImageSimpleModelTrackingScatterings(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        m.add_density_grid(np.array([[[1.e-30]]]), get_test_dust())

        s = m.add_point_source()
        s.name = 'first'
        s.luminosity = 1.
        s.temperature = 6000.

        s = m.add_point_source()
        s.name = 'second'
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_track_origin('scatterings', n_scat=5)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_image_invalid_option(self):

        # We can't use source_id and dust_id because tracking mode was not set
        # to 'detailed'

        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_image(source_id='all', component='source_emit')
        assert exc.value.args[0] == "cannot specify source_id since track_origin was not set to 'detailed'"

        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_image(dust_id='all', component='dust_emit')
        assert exc.value.args[0] == "cannot specify dust_id since track_origin was not set to 'detailed'"

        # The components should be 'source' and 'dust', anything else is invalid

        for component in ['source_emit', 'source_scat', 'dust_emit', 'dust_scat']:
            with pytest.raises(ValueError) as exc:
                wav, nufnu = self.m.get_image(n_scat=1, component=component)
            assert exc.value.args[0] == "component should be one of total/source/dust since track_origin='scatterings'"

    def test_image_n_scat_main_components(self):
        wav, nufnu = self.m.get_image(component='source')
        wav, nufnu = self.m.get_image(component='dust')

    def test_image_n_scat_n_scat_valid(self):
        for n_scat in range(6):
            wav, nufnu = self.m.get_image(n_scat=n_scat, component='source')
            wav, nufnu = self.m.get_image(n_scat=n_scat, component='dust')

    def test_image_n_scat_invalid(self):
        for n_scat in [-1, 6]:
            with pytest.raises(ValueError) as exc:
                wav, nufnu = self.m.get_image(n_scat=n_scat, component='source')
            assert exc.value.args[0] == 'n_scat should be between 0 and 5'

    def test_image_n_scat_values(self):
        for n_scat in range(6):
            image = self.m.get_image(n_scat=n_scat, component='source')
            if n_scat == 0:
                assert image.val.sum() > 0
            else:
                assert image.val.sum() == 0.


@pytest.mark.requires_hyperion_binaries
class TestSimpleModelInside(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_inside_observer((0., 0., 0.))
        i.set_image_limits(1., -1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_distance_fail(self):
        with pytest.raises(ValueError) as e:
            wav, nufnu = self.m.get_image(distance=1.)
        assert e.value.args[0] == 'Cannot specify distance for inside observers'


@pytest.mark.requires_hyperion_binaries
def test_regression_depth_bug(tmpdir):
    """
    This is a regression test for issue #21 reported by T. Bowers. If multiple
    images are requested with different depths, then if a photon did not fall
    in a depth interval, it was not included in subsequent image groups
    because 'return' was used instead 'cycle'.
    """

    m = Model()

    m.set_cartesian_grid([-1., 1.],
                         [-1., 1.],
                         [-1., 1.])

    m.add_density_grid(np.array([[[1.e-30]]]), get_test_dust())

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=False, image=True)
    i.set_viewing_angles([0.], [0.])
    i.set_image_limits(-1., 1., -1., 1.)
    i.set_image_size(1, 1)
    i.set_wavelength_range(1, 0.01, 1000.)
    i.set_depth(0.5, 1.0)

    i = m.add_peeled_images(sed=False, image=True)
    i.set_viewing_angles([0.], [0.])
    i.set_image_limits(-1., 1., -1., 1.)
    i.set_image_size(1, 1)
    i.set_wavelength_range(1, 0.01, 1000.)
    i.set_depth(-0.5, 0.5)

    i = m.add_peeled_images(sed=False, image=True)
    i.set_viewing_angles([0.], [0.])
    i.set_image_limits(-1., 1., -1., 1.)
    i.set_image_size(1, 1)
    i.set_wavelength_range(1, 0.01, 1000.)

    m.set_n_initial_iterations(0)

    m.set_n_photons(imaging=1)

    m.write(tmpdir.join(random_id()).strpath)

    mo = m.run(tmpdir.join(random_id()).strpath)

    wav, image1 = mo.get_image(group=0)
    wav, image2 = mo.get_image(group=1)
    wav, image3 = mo.get_image(group=2)

    assert image1.sum() == 0.
    assert image2.sum() > 0.
    assert image3.sum() == image2.sum()


@pytest.mark.requires_hyperion_binaries
class TestImage(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 2., -3., 4.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_depth(-2., 3.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_image_object(self):
        image = self.m.get_image(group=0)
        assert isinstance(image, Image)

    def test_image_attributes_no_distance(self):

        image = self.m.get_image(group=0, units='ergs/s')

        assert image.x_min == -1.
        assert image.x_max == +2.
        assert image.y_min == -3.
        assert image.y_max == +4.

        assert image.lon_min is None
        assert image.lon_max is None
        assert image.lat_min is None
        assert image.lat_max is None

        assert image.d_min == -2.
        assert image.d_max == 3.

        assert image.pix_area_sr is None

        assert image.distance is None

        assert not image.inside_observer

        assert image.units == 'ergs/s'

        assert image.nu.shape == (5,)
        assert image.wav.shape == (5,)
        assert image.val.shape == (2, 20, 10, 5)

    def test_image_attributes_distance(self):

        image = self.m.get_image(group=0, units='ergs/cm^2/s', distance=100.)

        assert image.x_min == -1.
        assert image.x_max == +2.
        assert image.y_min == -3.
        assert image.y_max == +4.

        lon_min = np.degrees(np.arctan(-1. / 100.))
        lon_max = np.degrees(np.arctan(+2. / 100.))
        lat_min = np.degrees(np.arctan(-3. / 100.))
        lat_max = np.degrees(np.arctan(+4. / 100.))

        assert_array_almost_equal_nulp(image.lon_min, lon_min, 5)
        assert_array_almost_equal_nulp(image.lon_max, lon_max, 5)
        assert_array_almost_equal_nulp(image.lat_min, lat_min, 5)
        assert_array_almost_equal_nulp(image.lat_max, lat_max, 5)

        pix_area_sr = np.radians(lon_max - lon_min) * np.radians(lat_max - lat_min) / 200

        assert_array_almost_equal_nulp(image.pix_area_sr, pix_area_sr, 5)

        assert image.distance == 100.

        assert not image.inside_observer

        assert image.units == 'ergs/cm^2/s'

        assert image.nu.shape == (5,)
        assert image.wav.shape == (5,)
        assert image.val.shape == (2, 20, 10, 5)

    def test_unit_conversion(self):

        # Assume that the initial scaling in ergs/cm^2/s is correct, so then
        # we just need to check the relative scaling.

        ref = self.m.get_image(group=0, units='ergs/cm^2/s', distance=100., inclination=1)

        # Make sure the flux is non-zero
        assert np.sum(ref.val) > 0

        # Check conversion to monochromatic flux
        mono = self.m.get_image(group=0, units='ergs/cm^2/s/Hz', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mono.val, 10)

        # Check conversion to Jy
        Jy = self.m.get_image(group=0, units='Jy', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), Jy.val * 1.e-23, 10)

        # Check conversion to mJy
        mJy = self.m.get_image(group=0, units='mJy', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mJy.val * 1.e-26, 10)

        # Check conversion to MJy/sr. For the far-field, all pixels have the
        # same area, so this is simple.
        MJy_per_sr = self.m.get_image(group=0, units='MJy/sr', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), MJy_per_sr.val * 1.e-17 * MJy_per_sr.pix_area_sr, 10)


@pytest.mark.requires_hyperion_binaries
class TestInsideImage(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        s = m.add_external_spherical_source()
        s.radius = 1.
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=False, image=True)
        i.set_inside_observer((0., 0., 0.))
        i.set_viewing_angles([1., 2., 3.], [1., 2., 3.])
        i.set_image_limits(5., -6., -7., 8.)
        i.set_image_size(30, 40)
        i.set_wavelength_range(3, 0.2, 50.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_image_object(self):
        image = self.m.get_image(group=0)
        assert isinstance(image, Image)

    def test_image_attributes_no_distance(self):

        with pytest.raises(ValueError) as exc:
            image = self.m.get_image(group=0, units='ergs/s')
        assert exc.value.args[0] == "Unknown units: ergs/s"

        image = self.m.get_image(group=0, units='ergs/cm^2/s')

        assert image.x_min is None
        assert image.x_max is None
        assert image.y_min is None
        assert image.y_max is None

        lon_min = +5.
        lon_max = -6.
        lat_min = -7.
        lat_max = +8.

        assert image.lon_min == lon_min
        assert image.lon_max == lon_max
        assert image.lat_min == lat_min
        assert image.lat_max == lat_max

        nx, ny = image.val.shape[2], image.val.shape[1]

        lon = np.linspace(np.radians(lon_min), np.radians(lon_max), nx + 1)
        lat = np.cos(np.linspace(np.radians(90. - lat_min), np.radians(90. - lat_max), ny + 1))
        dlon = lon[1:] - lon[:-1]
        dlat = lat[:-1] - lat[1:]
        DLON, DLAT = np.meshgrid(dlon, dlat)

        pix_area_sr = DLON * DLAT

        assert_array_almost_equal_nulp(image.pix_area_sr, pix_area_sr, 5)

        assert image.distance is None

        assert image.inside_observer

        assert image.units == 'ergs/cm^2/s'

        assert image.nu.shape == (3,)
        assert image.wav.shape == (3,)
        assert image.val.shape == (3, 40, 30, 3)

    def test_image_distance(self):

        with pytest.raises(ValueError) as exc:
            image = self.m.get_image(group=0, units='ergs/cm^2/s', distance=100.)
        assert exc.value.args[0] == "Cannot specify distance for inside observers"

    def test_unit_conversion(self):

        # Assume that the initial scaling in ergs/cm^2/s is correct, so then
        # we just need to check the relative scaling.

        ref = self.m.get_image(group=0, units='ergs/cm^2/s', inclination=0)

        # Make sure the flux is non-zero
        assert np.sum(ref.val) > 0

        # Check conversion to monochromatic flux
        mono = self.m.get_image(group=0, units='ergs/cm^2/s/Hz', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mono.val, 10)

        # Check conversion to Jy
        Jy = self.m.get_image(group=0, units='Jy', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), Jy.val * 1.e-23, 10)

        # Check conversion to mJy
        mJy = self.m.get_image(group=0, units='mJy', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mJy.val * 1.e-26, 10)

        # Check conversion to MJy/sr. For the far-field, all pixels have the
        # same area, so this is simple.
        MJy_per_sr = self.m.get_image(group=0, units='MJy/sr', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), MJy_per_sr.val * 1.e-17 * MJy_per_sr.pix_area_sr[:, :, np.newaxis], 10)


@pytest.mark.requires_hyperion_binaries
def test_flux_preserved_scatterings(tmpdir):
    # Regression test for issue #102 to ensure that flux is preserved when in
    # 'scatterings' mode

    dust = get_highly_reflective_dust()

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    m.add_density_grid(np.array([[[2.e-2]]]), dust)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 5000.

    i = m.add_peeled_images(sed=False, image=True)
    i.set_wavelength_range(10, 0.1, 1000.)
    i.set_viewing_angles([32], [33.])
    i.set_image_limits(-2., 2., -2., 2.)
    i.set_image_size(256, 256)
    i.set_track_origin('scatterings', n_scat=2)

    i = m.add_peeled_images(sed=False, image=True)
    i.set_wavelength_range(10, 0.1, 1000.)
    i.set_viewing_angles([32], [33.])
    i.set_image_limits(-2., 2., -2., 2.)
    i.set_image_size(256, 256)

    m.set_n_photons(initial=10000, imaging=10000)

    m.write(tmpdir.join('test.rtin').strpath)
    mo = m.run(tmpdir.join('test.rtout').strpath)

    image1 = mo.get_image(group=0, inclination=0)
    image2 = mo.get_image(group=1, inclination=0)

    np.testing.assert_allclose(image1.val, image2.val)

    source = mo.get_image(group=0, inclination=0, component='source')
    dust = mo.get_image(group=0, inclination=0, component='dust')

    np.testing.assert_allclose(image1.val, source.val + dust.val)


@pytest.mark.requires_hyperion_binaries
class TestImageStokesOption(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        img = m.add_peeled_images(sed=False, image=True)
        img.set_viewing_angles([1., 2.], [1., 2.])
        img.set_image_limits(-1., 2., -3., 4.)
        img.set_image_size(10, 20)
        img.set_wavelength_range(5, 0.1, 100.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()

        img.set_stokes(True)

        m.write(os.path.join(self.tmpdir, random_id()))
        self.m1 = m.run()

        img.set_stokes(False)

        m.write(os.path.join(self.tmpdir, random_id()))
        self.m2 = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_image_I(self):
        self.m1.get_image()
        self.m2.get_image()

    @pytest.mark.parametrize('stokes', ['Q', 'U', 'V', 'linpol', 'circpol'])
    def test_get_image_stokes(self, stokes):
        self.m1.get_image(stokes=stokes)
        with pytest.raises(ValueError) as exc:
            self.m2.get_image(stokes=stokes)
        assert exc.value.args[0] == "Only the Stokes I value was stored for this image"


def _inside_observer_direction(theta, phi):
    t, p = np.radians(theta), np.radians(phi)
    return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])


def _inside_observer_expected_lonlat(d, theta_v, phi_v):
    # Map (longitude, latitude) of a photon arriving from direction d, for an
    # observer looking towards (theta_v, phi_v). This is d expressed in the
    # local spherical basis (r_hat, phi_hat, -theta_hat) at the viewing
    # direction, converted to polar coordinates - i.e. an independent
    # reimplementation of the transform in images_peeled.f90.
    t, p = np.radians(theta_v), np.radians(phi_v)
    st, ct, cp, sp = np.sin(t), np.cos(t), np.cos(p), np.sin(p)
    vx = np.dot([st * cp, st * sp, ct], d)     # r_hat . d      -> map centre / +x
    vy = np.dot([-sp, cp, 0.], d)              # phi_hat . d    -> +longitude
    vz = np.dot([-ct * cp, -ct * sp, st], d)   # -theta_hat . d -> +latitude
    lon = np.degrees(np.arctan2(vy, vx))
    lat = np.degrees(np.arctan2(np.hypot(vx, vy), vz)) - 90.
    return lon, lat


@pytest.mark.requires_hyperion_binaries
def test_inside_observer_sky_coordinates(tmpdir):
    # Regression test for the inside-observer sky-coordinate transform. A point
    # source at cartesian position P, seen by an observer at the origin, sends
    # photons that arrive from direction d = -P/|P|. On the all-sky map for a
    # viewing direction (theta_v, phi_v) the source must appear at the position
    # obtained by expressing d in the local spherical frame of the viewing
    # direction. A previous implementation was discontinuous through
    # theta_v = 90 deg and ignored phi_v; the viewing angles below include
    # several cases it placed incorrectly (large phi at theta=90, either side
    # of theta=90, and general off-axis directions).

    # Place the source so that its photon arrival direction is +x.
    photon_dir = _inside_observer_direction(90., 0.)
    P = tuple(-0.5 * photon_dir)

    views = [(90., 0.), (90., 90.), (90., 150.), (60., 0.), (120., 0.),
             (89., 0.), (91., 0.), (45., 30.), (135., 200.)]

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    s = m.add_point_source()
    s.position = P
    s.luminosity = 1.
    s.temperature = 6000.
    i = m.add_peeled_images(sed=False, image=True)
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(180., -180., -90., 90.)
    i.set_image_size(360, 180)
    i.set_viewing_angles([v[0] for v in views], [v[1] for v in views])
    i.set_wavelength_range(1, 1., 1000.)
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=20000)
    m.set_seed(-1)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)

    val = out.get_image(group=0).val[:, :, :, 0]  # (n_view, n_y, n_x)
    n_y, n_x = val.shape[1], val.shape[2]
    xmin, xmax, ymin, ymax = 180., -180., -90., 90.

    for iv, (theta_v, phi_v) in enumerate(views):
        img = val[iv]
        ys, xs = np.nonzero(img > 0)
        assert len(xs) > 0, "source not found for view (%.0f, %.0f)" % (theta_v, phi_v)
        w = img[ys, xs]
        lon = ((xmin + (xs + 0.5) * (xmax - xmin) / n_x) * w).sum() / w.sum()
        lat = ((ymin + (ys + 0.5) * (ymax - ymin) / n_y) * w).sum() / w.sum()

        # (1) Convention-independent invariant: the great-circle distance on the
        # map from the centre (0, 0) to the source must equal the true angle
        # between the photon direction and the viewing direction (any rigid
        # rotation preserves angles). This catches the tearing and the ignored
        # phi of the old code.
        map_sep = np.degrees(np.arccos(np.clip(
            np.cos(np.radians(lat)) * np.cos(np.radians(lon)), -1., 1.)))
        true_sep = np.degrees(np.arccos(np.clip(
            np.dot(photon_dir, _inside_observer_direction(theta_v, phi_v)), -1., 1.)))
        assert abs(map_sep - true_sep) < 2., \
            "view (%.0f, %.0f): map separation %.1f deg != true %.1f deg" % (
                theta_v, phi_v, map_sep, true_sep)

        # (2) Full position matches the independent reference (catches roll and
        # reflection errors that preserve the angular separation).
        exp_lon, exp_lat = _inside_observer_expected_lonlat(photon_dir, theta_v, phi_v)
        dlon = (lon - exp_lon + 180.) % 360. - 180.
        assert abs(dlon) < 2. and abs(lat - exp_lat) < 2., \
            "view (%.0f, %.0f): (%.1f, %.1f) != expected (%.1f, %.1f)" % (
                theta_v, phi_v, lon, lat, exp_lon, exp_lat)


def test_image_and_sed_units_constructor():
    # Regression test: the ``units`` constructor argument of Image and SED must
    # be routed through the validated ``unit`` property, so that it is stored
    # (readable via ``.unit``) and invalid values are rejected. Previously it
    # was assigned to an unvalidated attribute and left ``.unit`` unset.
    from ..sed import SED
    for cls in (Image, SED):
        obj = cls(nu=np.array([1.e10, 1.e11]), units='Jy')
        assert obj.unit == 'Jy'
        assert obj.units == 'Jy'
        with pytest.raises(ValueError):
            cls(nu=np.array([1.e10, 1.e11]), units=42)


@pytest.mark.requires_hyperion_binaries
def test_inside_observer_flux_dilution(tmpdir):
    # Regression test for the inside-observer 1/d^2 flux dilution. Two point
    # sources of equal luminosity at distances d1 and d2 from the observer must
    # have observed fluxes in the ratio (d2/d1)**2. A previous implementation
    # divided by (d - d_min)**2 instead of d**2, giving the wrong ratio whenever
    # a nonzero depth minimum was set.
    d1, d2, d_min = 2., 5., 1.
    m = Model()
    m.set_cartesian_grid([-8., 8.], [-8., 8.], [-8., 8.])
    for pos in [(d1, 0., 0.), (0., d2, 0.)]:
        s = m.add_point_source()
        s.position = pos
        s.luminosity = 1.
        s.temperature = 6000.
    i = m.add_peeled_images(sed=False, image=True)
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(180., -180., -90., 90.)
    i.set_image_size(360, 180)
    i.set_viewing_angles([90.], [0.])
    i.set_wavelength_range(1, 1., 1000.)
    i.set_depth(d_min, 20.)   # nonzero d_min exposes the bug
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=200000)
    m.set_seed(-1)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)

    val = out.get_image(group=0).val[0, :, :, 0]   # single view, single wavelength
    brightest = np.sort(val[val > 0])[::-1]
    assert len(brightest) >= 2, "expected two point sources in the image"
    # the closer source (d1) is the brighter one
    ratio = brightest[0] / brightest[1]
    np.testing.assert_allclose(ratio, (d2 / d1) ** 2, rtol=0.08)


@pytest.mark.requires_hyperion_binaries
def test_inside_observer_peeloff_optical_depth(tmpdir):
    # Regression test for integrating the peeloff optical depth all the way to
    # the (inside) observer. A point source at distance d sits in uniform,
    # purely-absorbing dust of flat opacity chi and density rho, so its direct
    # light is attenuated by exp(-chi*rho*d) over the full path to the observer.
    # A previous implementation integrated the optical depth only to d - d_min,
    # under-attenuating the source whenever a nonzero depth minimum was set
    # (comparing the dust/no-dust flux ratio at fixed distance cancels the
    # 1/d^2 dilution, isolating the optical-depth path).
    from ...dust import IsotropicDust

    d, chi, rho, d_min = 1.e16, 1., 1.e-16, 0.5e16
    tau_full = chi * rho * d   # = 1.0

    def peak_flux(with_dust):
        m = Model()
        m.set_cartesian_grid([-2.e16, 2.e16], [-2.e16, 2.e16], [-2.e16, 2.e16])
        if with_dust:
            dust = IsotropicDust([3.e9, 3.e16], [0., 0.], [chi, chi])
            dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)
            m.add_density_grid(np.ones((1, 1, 1)) * rho, dust)
        s = m.add_point_source()
        s.position = (d, 0., 0.)
        s.luminosity = 1.
        s.temperature = 6000.
        i = m.add_peeled_images(sed=False, image=True)
        i.set_inside_observer((0., 0., 0.))
        i.set_image_limits(180., -180., -90., 90.)
        i.set_image_size(360, 180)
        i.set_viewing_angles([90.], [0.])
        i.set_wavelength_range(1, 1., 1000.)
        i.set_depth(d_min, 3.e16)   # nonzero d_min previously truncated the tau path
        m.set_n_initial_iterations(0)
        m.set_n_photons(imaging=100000)
        m.set_seed(-1)
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        return out.get_image(group=0).val[0, :, :, 0].max()

    ratio = peak_flux(with_dust=True) / peak_flux(with_dust=False)
    np.testing.assert_allclose(ratio, np.exp(-tau_full), rtol=0.05)
