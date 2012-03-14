import numpy as np
from hyperion.model import Model
from hyperion.util.functions import random_filename

import tempfile
import os
import string
import random
import pytest


class TestSimpleModel(object):

    def setup_class(self):

        m = Model()

        m.set_cartesian_grid([-1., 1.],
                             [-1., 1.],
                             [-1., 1.])

        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 6000.

        i = m.add_peeled_images()
        i.set_viewing_angles([1., 2.], [1., 2.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_image_size(10, 20)
        i.set_wavelength_range(5, 0.1, 100.)
        i.set_aperture_range(3, 1., 10.)

        m.set_n_initial_iterations(0)

        m.set_n_photons(imaging=1)

        m.write(random_filename())

        self.m = m.run()

    @pytest.mark.xfail
    def test_sed_group(self):
        wav, nufnu = self.m.get_sed(group=0)

    def test_sed_group_invalid1(self):
        with pytest.raises(Exception):
            wav, nufnu = self.m.get_sed(group=-2)
            # negative indexing allowed, but only one group present

    @pytest.mark.xfail
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
