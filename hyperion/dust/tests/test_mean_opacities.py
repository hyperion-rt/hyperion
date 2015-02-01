from __future__ import print_function, division

import matplotlib.pyplot as plt

import numpy as np
from numpy.testing import assert_allclose
from astropy.tests.helper import pytest

from ..optical_properties import OpticalProperties
from ..mean_opacities import MeanOpacities
from ...util.functions import virtual_file, B_nu
from ...util.constants import sigma


class TestMeanOpacities(object):

    def setup_class(self):
        self.o = OpticalProperties()
        self.o.nu = np.array([1.e8, 1.e16])
        self.o.chi = np.array([1.e-2, 1])
        self.o.albedo = np.array([0., 0.5])

    def test_init(self):
        MeanOpacities()

    def test_compute(self):

        m = MeanOpacities()
        m.compute(self.o, n_temp=10, temp_min=1., temp_max=1000.)

        assert_allclose(m.temperature[0], 1.)
        assert_allclose(m.temperature[-1], 1000.)

        assert_allclose(m.specific_energy,
                                   4. * sigma * m.temperature ** 4. * m.kappa_planck)

    def test_io(self):

        m = MeanOpacities()
        m.compute(self.o, n_temp=10, temp_min=1., temp_max=1000.)

        f = virtual_file()
        m.to_hdf5_group(f)
        m_new = MeanOpacities()
        m_new.from_hdf5_group(f)

        assert_allclose(m.temperature, m_new.temperature)
        assert_allclose(m.specific_energy, m_new.specific_energy)
        assert_allclose(m.chi_planck, m_new.chi_planck)
        assert_allclose(m.kappa_planck, m_new.kappa_planck)
        assert_allclose(m.chi_inv_planck, m_new.chi_inv_planck)
        assert_allclose(m.kappa_inv_planck, m_new.kappa_inv_planck)
        assert_allclose(m.chi_rosseland, m_new.chi_rosseland)
        assert_allclose(m.kappa_rosseland, m_new.kappa_rosseland)
        assert m.hash() == m_new.hash()

    def test_plot(self):

        # Just check that plot runs without crashing

        fig = plt.figure()

        m = MeanOpacities()
        m.compute(self.o, n_temp=10, temp_min=1., temp_max=1000.)
        m.plot(fig, 111)

        plt.close(fig)
