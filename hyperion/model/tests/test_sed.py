from __future__ import print_function, division

import pytest
import numpy as np

from .. import Model
from ...util.functions import random_filename
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

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        m.write(random_filename())

        self.m = m.run()

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
        s.luminosity = 1.
        s.temperature = 6000.

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images(sed=True, image=False)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)
        i.set_track_origin('detailed')

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        m.write(random_filename())

        self.m = m.run()

    def test_sed_source_all(self):
        wav, nufnu = self.m.get_sed(source_id='all', component='source_emit')

    def test_sed_source_valid1(self):
        wav, nufnu = self.m.get_sed(source_id=0, component='source_emit')

    def test_sed_source_valid2(self):
        wav, nufnu = self.m.get_sed(source_id=1, component='source_emit')

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

        m.write(random_filename())

        self.m = m.run()

    def test_distance_fail(self):
        with pytest.raises(ValueError) as e:
            wav, nufnu = self.m.get_sed(distance=1.)
        assert e.value.args[0] == 'Cannot specify distance for inside observers'
