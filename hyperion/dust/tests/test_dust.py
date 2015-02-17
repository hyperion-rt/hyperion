import numpy as np
from numpy.testing import assert_allclose
from astropy.tests.helper import pytest
import matplotlib.pyplot as plt

from .. import SphericalDust, IsotropicDust
from ...util.functions import random_id, B_nu


def test_missing_properties(tmpdir):
    d = SphericalDust()
    with pytest.raises(Exception) as e:
        d.write(tmpdir.join(random_id()).strpath)
    assert e.value.args[0] == "The following attributes of the optical properties have not been set: nu, chi, albedo, mu, P1, P2, P3, P4"


class TestSphericalDust(object):

    def setup_method(self, method):

        self.dust = SphericalDust()

        self.dust.optical_properties.nu = np.logspace(0., 20., 100)
        self.dust.optical_properties.albedo = np.repeat(0.5, 100)
        self.dust.optical_properties.chi = np.ones(100)
        self.dust.optical_properties.mu = [-1., 1.]
        self.dust.optical_properties.initialize_scattering_matrix()

    def test_helpers(self):

        nu = np.logspace(5., 15., 1000)

        # Here we don't set the mean opacities to make sure they are computed
        # automatically

        assert_allclose(self.dust.kappa_nu_temperature(34.),
                        self.dust.kappa_nu_spectrum(nu, B_nu(nu, 34)))

        assert_allclose(self.dust.chi_nu_temperature(34.),
                        self.dust.chi_nu_spectrum(nu, B_nu(nu, 34)))

    def test_conversions_out_of_bounds(self):

        np.random.seed(12345)

        self.dust.mean_opacities.compute(self.dust.optical_properties,
                                         n_temp=10, temp_min=1., temp_max=1000.)

        # Test with scalars
        assert_allclose(self.dust.temperature2specific_energy(0.1),
                                   self.dust.mean_opacities.specific_energy[0])
        assert_allclose(self.dust.temperature2specific_energy(1e4),
                                   self.dust.mean_opacities.specific_energy[-1])

        # Test with arrays
        assert_allclose(self.dust.temperature2specific_energy(np.array([0.1])),
                                   self.dust.mean_opacities.specific_energy[0])
        assert_allclose(self.dust.temperature2specific_energy(np.array([1e4])),
                                   self.dust.mean_opacities.specific_energy[-1])

        # Test with scalars
        assert_allclose(self.dust.specific_energy2temperature(1.e-10),
                                   self.dust.mean_opacities.temperature[0])
        assert_allclose(self.dust.specific_energy2temperature(1.e+10),
                                   self.dust.mean_opacities.temperature[-1])

        # Test with arrays
        assert_allclose(self.dust.specific_energy2temperature(np.array([1.e-10])),
                                   self.dust.mean_opacities.temperature[0])
        assert_allclose(self.dust.specific_energy2temperature(np.array([1.e+10])),
                                   self.dust.mean_opacities.temperature[-1])

    def test_conversions_roundtrip(self):

        np.random.seed(12345)

        # Here we don't set the mean opacities to make sure they are computed
        # automatically

        temperatures1 = 10. ** np.random.uniform(0., 3., 100)

        specific_energies = self.dust.temperature2specific_energy(temperatures1)

        temperatures2 = self.dust.specific_energy2temperature(specific_energies)

        assert_allclose(temperatures1, temperatures2)

    def test_set_lte(self):

        np.random.seed(12345)

        # Here we don't set the mean opacities to make sure they are computed
        # automatically

        self.dust.set_lte_emissivities(n_temp=10, temp_min=1., temp_max=1000.)

        assert_allclose(self.dust.mean_opacities.temperature[0], 1)
        assert_allclose(self.dust.mean_opacities.temperature[-1], 1000)

        assert_allclose(self.dust.mean_opacities.specific_energy,
                                  self.dust.emissivities.var)

    def test_plot(self, tmpdir):

        # Here we don't set the mean opacities or the emissivities to make sure
        # they are computed automatically

        self.dust.plot(tmpdir.join('test.png').strpath)

    def test_hash(self, tmpdir):

        # Here we don't set the mean opacities or the emissivities to make sure
        # they are computed automatically

        assert self.dust.hash() == '8e1f63eedcafcc05183b99ee8fe1333a'

    def test_io(self, tmpdir):

        filename = tmpdir.join('test.hdf5').strpath

        self.dust.write(filename)

        dust2 = SphericalDust(filename)

        assert_allclose(self.dust.optical_properties.nu, dust2.optical_properties.nu)
        assert_allclose(self.dust.optical_properties.chi, dust2.optical_properties.chi)
        assert_allclose(self.dust.optical_properties.albedo, dust2.optical_properties.albedo)
        assert_allclose(self.dust.optical_properties.mu, dust2.optical_properties.mu)
        assert_allclose(self.dust.optical_properties.P1, dust2.optical_properties.P1)
        assert_allclose(self.dust.optical_properties.P2, dust2.optical_properties.P2)
        assert_allclose(self.dust.optical_properties.P3, dust2.optical_properties.P3)
        assert_allclose(self.dust.optical_properties.P4, dust2.optical_properties.P4)

        assert_allclose(self.dust.mean_opacities.temperature, dust2.mean_opacities.temperature)
        assert_allclose(self.dust.mean_opacities.specific_energy, dust2.mean_opacities.specific_energy)
        assert_allclose(self.dust.mean_opacities.chi_planck, dust2.mean_opacities.chi_planck)
        assert_allclose(self.dust.mean_opacities.kappa_planck, dust2.mean_opacities.kappa_planck)
        assert_allclose(self.dust.mean_opacities.chi_inv_planck, dust2.mean_opacities.chi_inv_planck)
        assert_allclose(self.dust.mean_opacities.kappa_inv_planck, dust2.mean_opacities.kappa_inv_planck)
        assert_allclose(self.dust.mean_opacities.chi_rosseland, dust2.mean_opacities.chi_rosseland)
        assert_allclose(self.dust.mean_opacities.kappa_rosseland, dust2.mean_opacities.kappa_rosseland)

        assert self.dust.emissivities.is_lte == dust2.emissivities.is_lte
        assert self.dust.emissivities.var_name == dust2.emissivities.var_name
        assert_allclose(self.dust.emissivities.nu, dust2.emissivities.nu)
        assert_allclose(self.dust.emissivities.var, dust2.emissivities.var)
        assert_allclose(self.dust.emissivities.jnu, dust2.emissivities.jnu)

        assert self.dust.sublimation_mode == dust2.sublimation_mode
        assert_allclose(self.dust.sublimation_energy, dust2.sublimation_energy)


def test_isotropic_dust():
    nu = np.logspace(0., 20., 100)
    albedo = np.repeat(0.5, 100)
    chi = np.ones(100)
    dust = IsotropicDust(nu, albedo, chi)
