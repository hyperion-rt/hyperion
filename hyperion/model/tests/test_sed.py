from __future__ import print_function, division

import os
import shutil
import tempfile

from astropy.tests.helper import pytest
import numpy as np

from numpy.testing import assert_array_almost_equal_nulp

from .. import Model
from ..sed import SED
from ...util.functions import random_id
from .test_helpers import get_test_dust


class TestSEDSimpleModel(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=True, image=False)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)
        i.set_stokes(True)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_sed_group(self):
        wav, nufnu = self.m.get_sed(group=0)

    def test_sed_group_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(group=-2)
            # negative indexing allowed, but only one group present
        assert exc.value.args[0] == 'File only contains 1 image/SED group(s)'

    def test_sed_group_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(group=1)
            # zero-based, and only one group present
        assert exc.value.args[0] == 'File only contains 1 image/SED group(s)'

    def test_sed_dim(self):
        wav, nufnu = self.m.get_sed()
        assert nufnu.shape == (2, 3, 5)

    def test_sed_dim_incl1(self):
        wav, nufnu = self.m.get_sed(inclination=0)
        assert nufnu.shape == (3, 5)

    def test_sed_dim_incl2(self):
        wav, nufnu = self.m.get_sed(inclination=1)
        assert nufnu.shape == (3, 5)

    def test_sed_dim_incl_invalid1(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_sed(inclination=2)

    def test_sed_dim_incl_invalid2(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_sed(inclination=-3)

    def test_sed_dim_incl_invalid3(self):
        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_sed(inclination=12.3)
        assert exc.value.args[0] == "inclination should be an integer (it should be the index of the inclination, not the value itself)"

    def test_sed_dim_aper1(self):
        wav, nufnu = self.m.get_sed(aperture=0)
        assert nufnu.shape == (2, 5)

    def test_sed_dim_aper2(self):
        wav, nufnu = self.m.get_sed(aperture=2)
        assert nufnu.shape == (2, 5)

    def test_sed_dim_aper_invalid1(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_sed(aperture=3)

    def test_sed_dim_aper_invalid2(self):
        with pytest.raises(IndexError):
            wav, nufnu = self.m.get_sed(aperture=-4)

    def test_sed_dim_aper_invalid3(self):
        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_sed(aperture=344.3)
        assert exc.value.args[0] == "aperture should be an integer (it should be the index of the aperture, not the value itself)"

    @pytest.mark.parametrize(('stokes'), ['I', 'Q', 'U', 'V',
                                          'linpol', 'circpol'])
    def test_sed_stokes(self, stokes):
        wav, nufnu = self.m.get_sed(stokes=stokes)
        assert nufnu.shape == (2, 3, 5)

    @pytest.mark.parametrize(('stokes'), ['A', 'b', 1, (3,),  # invalid values
                                          'i', 'q', 'u', 'v'])  # lowercase
    def test_sed_stokes_invalid(self, stokes):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(stokes=stokes)
        if isinstance(stokes, basestring):
            assert exc.value.args[0] == "Unknown Stokes parameter: %s" % stokes
        else:
            assert exc.value.args[0] == "stokes argument should be a string"

    @pytest.mark.parametrize(('units'), ['ergs/s'])
    def test_sed_nodistance_units(self, units):
        wav, nufnu = self.m.get_sed(units=units)

    @pytest.mark.parametrize(('units'), ['ergs/cm^2/s', 'mJy', 'Jy', 'ergs/cm^2/s/Hz', 'MJy/sr'])
    def test_sed_nodistance_units_invalid(self, units):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(units=units)
        assert exc.value.args[0] == 'Since distance= is not specified, units should be set to ergs/s'


class TestSEDSimpleModelTrackingDetailed(object):

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

        i = m.add_peeled_images(sed=True, image=False)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)
        i.set_track_origin('detailed')

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_sed_source_all(self):
        wav, nufnu = self.m.get_sed(source_id='all', component='source_emit')

    def test_sed_source_valid1(self):
        wav, nufnu = self.m.get_sed(source_id=0, component='source_emit')

    def test_sed_source_valid2(self):
        wav, nufnu = self.m.get_sed(source_id=1, component='source_emit')

    def test_sed_source_valid3(self):
        wav, nufnu = self.m.get_sed(source_id='first', component='source_emit')

    def test_sed_source_valid4(self):
        wav, nufnu = self.m.get_sed(source_id='second', component='source_emit')

    def test_sed_source_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(source_id=-1, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

    def test_sed_source_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(source_id=2, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

    def test_sed_dust_all(self):
        wav, nufnu = self.m.get_sed(dust_id='all', component='dust_emit')

    def test_sed_dust_valid1(self):
        wav, nufnu = self.m.get_sed(dust_id=0, component='dust_emit')

    def test_sed_dust_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(dust_id=-1, component='dust_emit')
        assert exc.value.args[0] == 'dust_id should be between 0 and 0'

    def test_sed_dust_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_sed(dust_id=1, component='dust_emit')
        assert exc.value.args[0] == 'dust_id should be between 0 and 0'


class TestSEDSimpleModelTrackingScatterings(object):

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

        i = m.add_peeled_images(sed=True, image=False)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)
        i.set_track_origin('scatterings', n_scat=5)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_sed_invalid_option(self):

        # We can't use source_id and dust_id because tracking mode was not set
        # to 'detailed'

        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_sed(source_id='all', component='source_emit')
        assert exc.value.args[0] == "cannot specify source_id since track_origin was not set to 'detailed'"

        with pytest.raises(Exception) as exc:
            wav, nufnu = self.m.get_sed(dust_id='all', component='dust_emit')
        assert exc.value.args[0] == "cannot specify dust_id since track_origin was not set to 'detailed'"

        # The components should be 'source' and 'dust', anything else is invalid

        for component in ['source_emit', 'source_scat', 'dust_emit', 'dust_scat']:
            with pytest.raises(ValueError) as exc:
                wav, nufnu = self.m.get_sed(n_scat=1, component=component)
            assert exc.value.args[0] == "component should be one of total/source/dust since track_origin='scatterings'"

    def test_sed_n_scat_main_components(self):
        wav, nufnu = self.m.get_sed(component='source')
        wav, nufnu = self.m.get_sed(component='dust')

    def test_sed_n_scat_n_scat_valid(self):
        for n_scat in range(6):
            wav, nufnu = self.m.get_sed(n_scat=n_scat, component='source')
            wav, nufnu = self.m.get_sed(n_scat=n_scat, component='dust')

    def test_sed_n_scat_invalid(self):
        for n_scat in [-1, 6]:
            with pytest.raises(ValueError) as exc:
                wav, nufnu = self.m.get_sed(n_scat=n_scat, component='source')
            assert exc.value.args[0] == 'n_scat should be between 0 and 5'

    def test_sed_n_scat_values(self):
        for n_scat in range(6):
            sed = self.m.get_sed(n_scat=n_scat, component='source')
            if n_scat == 0:
                assert sed.val.sum() > 0
            else:
                assert sed.val.sum() == 0.

class TestSimpleModelInside(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=True, image=False)
        i.set_inside_observer((0., 0., 0.))
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_distance_fail(self):
        with pytest.raises(ValueError) as e:
            wav, nufnu = self.m.get_sed(distance=1.)
        assert e.value.args[0] == 'Cannot specify distance for inside observers'


class TestSED(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        sed = m.add_peeled_images(sed=True, image=False)
        sed.set_viewing_angles([1., 2.], [1., 2.])
        sed.set_wavelength_range(5, 0.1, 100.)
        sed.set_aperture_range(4, 2., 5.)
        sed.set_depth(-2., 3.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_sed_object(self):
        sed = self.m.get_sed(group=0)
        assert isinstance(sed, SED)

    def test_sed_attributes_no_distance(self):

        sed = self.m.get_sed(group=0, units='ergs/s')

        assert sed.ap_min == 2.
        assert sed.ap_max == 5.

        assert sed.d_min == -2.
        assert sed.d_max == 3.

        assert sed.distance is None

        assert not sed.inside_observer

        assert sed.units == 'ergs/s'

        assert sed.nu.shape == (5,)
        assert sed.wav.shape == (5,)
        assert sed.val.shape == (2, 4, 5)

    def test_sed_attributes_distance(self):

        sed = self.m.get_sed(group=0, units='ergs/cm^2/s', distance=100.)

        assert sed.ap_min == 2.
        assert sed.ap_max == 5.

        assert sed.distance == 100.

        assert not sed.inside_observer

        assert sed.units == 'ergs/cm^2/s'

        assert sed.nu.shape == (5,)
        assert sed.wav.shape == (5,)
        assert sed.val.shape == (2, 4, 5)

    def test_unit_conversion(self):

        # Assume that the initial scaling in ergs/cm^2/s is correct, so then
        # we just need to check the relative scaling.

        ref = self.m.get_sed(group=0, units='ergs/cm^2/s', distance=100., inclination=1)

        # Make sure the flux is non-zero
        assert np.sum(ref.val) > 0

        # Check conversion to monochromatic flux
        mono = self.m.get_sed(group=0, units='ergs/cm^2/s/Hz', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mono.val, 10)

        # Check conversion to Jy
        Jy = self.m.get_sed(group=0, units='Jy', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), Jy.val * 1.e-23, 10)

        # Check conversion to mJy
        mJy = self.m.get_sed(group=0, units='mJy', distance=100., inclination=1)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mJy.val * 1.e-26, 10)


class TestInsideSED(object):

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

        sed = m.add_peeled_images(sed=True, image=False)
        sed.set_inside_observer((0., 0., 0.))
        sed.set_viewing_angles([1., 2., 3.], [1., 2., 3.])
        sed.set_wavelength_range(3, 0.2, 50.)
        sed.set_aperture_range(4, 2., 50.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()
        m.write(os.path.join(self.tmpdir, random_id()))

        self.m = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_sed_object(self):
        sed = self.m.get_sed(group=0)
        assert isinstance(sed, SED)

    def test_sed_attributes_no_distance(self):

        with pytest.raises(ValueError) as exc:
            sed = self.m.get_sed(group=0, units='ergs/s')
        assert exc.value.args[0] == "Unknown units: ergs/s"

        sed = self.m.get_sed(group=0, units='ergs/cm^2/s')

        assert sed.ap_min == 2.
        assert sed.ap_max == 50.

        assert sed.distance is None

        assert sed.inside_observer

        assert sed.units == 'ergs/cm^2/s'

        assert sed.nu.shape == (3,)
        assert sed.wav.shape == (3,)
        assert sed.val.shape == (3, 4, 3)

    def test_sed_distance(self):

        with pytest.raises(ValueError) as exc:
            sed = self.m.get_sed(group=0, units='ergs/cm^2/s', distance=100.)
        assert exc.value.args[0] == "Cannot specify distance for inside observers"

    def test_unit_conversion(self):

        # Assume that the initial scaling in ergs/cm^2/s is correct, so then
        # we just need to check the relative scaling.

        ref = self.m.get_sed(group=0, units='ergs/cm^2/s', inclination=0)

        # Make sure the flux is non-zero
        assert np.sum(ref.val) > 0

        # Check conversion to monochromatic flux
        mono = self.m.get_sed(group=0, units='ergs/cm^2/s/Hz', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mono.val, 10)

        # Check conversion to Jy
        Jy = self.m.get_sed(group=0, units='Jy', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), Jy.val * 1.e-23, 10)

        # Check conversion to mJy
        mJy = self.m.get_sed(group=0, units='mJy', inclination=0)
        assert_array_almost_equal_nulp((ref.val / ref.nu), mJy.val * 1.e-26, 10)


class TestSEDStokesOption(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        sed = m.add_peeled_images(sed=True, image=False)
        sed.set_viewing_angles([1., 2.], [1., 2.])
        sed.set_wavelength_range(5, 0.1, 100.)
        sed.set_aperture_range(4, 2., 5.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=10000)

        self.tmpdir = tempfile.mkdtemp()

        sed.set_stokes(True)

        m.write(os.path.join(self.tmpdir, random_id()))
        self.m1 = m.run()

        sed.set_stokes(False)

        m.write(os.path.join(self.tmpdir, random_id()))
        self.m2 = m.run()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_get_sed_I(self):
        self.m1.get_sed()
        self.m2.get_sed()

    @pytest.mark.parametrize('stokes', ['Q', 'U', 'V', 'linpol', 'circpol'])
    def test_get_sed_stokes(self, stokes):
        self.m1.get_sed(stokes=stokes)
        with pytest.raises(ValueError) as exc:
            self.m2.get_sed(stokes=stokes)
        assert exc.value.args[0] == "Only the Stokes I value was stored for this SED"
