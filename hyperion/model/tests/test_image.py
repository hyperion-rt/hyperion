from __future__ import print_function, division

import pytest
import numpy as np

from .. import Model
from ...util.functions import random_filename
from .test_helpers import get_test_dust


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

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        m.write(random_filename())

        self.m = m.run()

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

        i = m.add_peeled_images(sed=False, image=True)
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_track_origin('detailed')

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        m.write(random_filename())

        self.m = m.run()

    def test_image_source_all(self):
        wav, nufnu = self.m.get_image(source_id='all', component='source_emit')

    def test_image_source_valid1(self):
        wav, nufnu = self.m.get_image(source_id=0, component='source_emit')

    def test_image_source_valid2(self):
        wav, nufnu = self.m.get_image(source_id=1, component='source_emit')

    def test_image_source_invalid1(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(source_id=-1, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

    def test_image_source_invalid2(self):
        with pytest.raises(ValueError) as exc:
            wav, nufnu = self.m.get_image(source_id=2, component='source_emit')
        assert exc.value.args[0] == 'source_id should be between 0 and 1'

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

        m.write(random_filename())

        self.m = m.run()

    def test_distance_fail(self):
        with pytest.raises(ValueError) as e:
            wav, nufnu = self.m.get_image(distance=1.)
        assert e.value.args[0] == 'Cannot specify distance for inside observers'


def test_regression_depth_bug():
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

    m.write(random_filename())

    mo = m.run(random_filename())

    wav, image1 = mo.get_image(group=0)
    wav, image2 = mo.get_image(group=1)
    wav, image3 = mo.get_image(group=2)

    assert image1.sum() == 0.
    assert image2.sum() > 0.
    assert image3.sum() == image2.sum()
