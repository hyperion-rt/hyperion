from __future__ import print_function, division

import os
import random
import string
import tempfile

import pytest
import numpy as np

from .. import Model
from ...util.functions import random_filename


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
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(group=-2)
            # negative indexing allowed, but only one group present

    def test_sed_group_invalid2(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(group=1)
            # zero-based, and only one group present

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
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(inclination=2)

    def test_sed_dim_incl_invalid2(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(inclination=-3)

    def test_sed_dim_aper1(self):
        wav, nufnu = self.m.get_sed(aperture=0)
        assert nufnu.shape == (2, 5)

    def test_sed_dim_aper2(self):
        wav, nufnu = self.m.get_sed(aperture=2)
        assert nufnu.shape == (2, 5)

    def test_sed_dim_aper_invalid1(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(aperture=3)

    def test_sed_dim_aper_invalid2(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(aperture=-4)

    @pytest.mark.parametrize(('stokes'), ['I', 'Q', 'U', 'V',
                                          'linpol', 'circpol'])
    def test_sed_stokes(self, stokes):
        wav, nufnu = self.m.get_sed(stokes=stokes)
        assert nufnu.shape == (2, 3, 5)

    @pytest.mark.parametrize(('stokes'), ['A', 'b', 1, (3,),  # invalid values
                                          'i', 'q', 'u', 'v'])  # lowercase
    def test_sed_stokes_invalid(self, stokes):
        with pytest.raises(ValueError):
            wav, nufnu = self.m.get_sed(stokes=stokes)

    @pytest.mark.parametrize(('units'), ['ergs/s'])
    def test_sed_nodistance_units(self, units):
        wav, nufnu = self.m.get_sed(units=units)

    @pytest.mark.parametrize(('units'), ['ergs/cm^2/s', 'mJy', 'Jy', 'ergs/cm^2/s/Hz', 'MJy/sr'])
    def test_sed_nodistance_units_invalid(self, units):
        with pytest.raises(ValueError):
            wav, nufnu = self.m.get_sed(units=units)


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
