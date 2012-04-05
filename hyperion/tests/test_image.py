import tempfile
import os
import string
import random

import numpy as np
import pytest

from hyperion.model import Model
from hyperion.util.functions import random_filename


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

    @pytest.mark.xfail
    def test_image_group(self):
        wav, nufnu = self.m.get_image(group=0)

    def test_image_group_invalid1(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_image(group=-2)
            # negative indexing allowed, but only one group present

    @pytest.mark.xfail
    def test_image_group_invalid2(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_image(group=1)
            # zero-based, and only one group present

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
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_image(inclination=2)

    def test_image_dim_incl_invalid2(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_image(inclination=-3)

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
        i.set_inside_observer((0.,0.,0.))
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
        assert e.value.message == 'Cannot specify distance for inside observers'
