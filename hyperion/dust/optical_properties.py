from __future__ import print_function, division

import atpy
import numpy as np

from ..util.integrate import integrate_loglog, integrate_linlog_subset
from ..util.interpolate import interp1d_fast, interp1d_fast_loglog, \
                                      interp1d_fast_linlog
from ..util.functions import extrap1d_log10, B_nu, FreezableClass, \
                                    nu_common, planck_nu_range, is_numpy_array, monotonically_increasing
from ..util.constants import c, sigma
from ..util.logger import logger


class OpticalProperties(FreezableClass):

    def __init__(self):

        # Wavelengths
        self.nu = None

        # Opacity to extinction
        self.chi = None

        # Albedo
        self.albedo = None

        # Scattering angles
        self.mu = None

        # Scattering matrix elements
        self.P1 = None
        self.P2 = None
        self.P3 = None
        self.P4 = None

        # Prevent new attributes
        self._freeze()

    def __getattr__(self, attribute):
        if attribute == 'kappa':
            return self.chi * (1. - self.albedo)
        elif attribute == 'wav':
            return c / self.nu * 1.e4
        else:
            raise AttributeError(attribute)

    def _sort(self):

        if self.mu[-1] < self.mu[0]:

            self.mu = self.mu[::-1]
            self.P1 = self.P1[:, ::-1]
            self.P2 = self.P2[:, ::-1]
            self.P3 = self.P3[:, ::-1]
            self.P4 = self.P4[:, ::-1]

        if self.nu[-1] < self.nu[0]:

            self.wav = self.wav[::-1]
            self.nu = self.nu[::-1]
            self.albedo = self.albedo[::-1]
            self.chi = self.chi[::-1]
            self.P1 = self.P1[::-1, :]
            self.P2 = self.P2[::-1, :]
            self.P3 = self.P3[::-1, :]
            self.P4 = self.P4[::-1, :]

    def initialize_scattering_matrix(self):

        self.P1 = np.zeros((len(self.wav), len(self.mu)))
        self.P2 = np.zeros((len(self.wav), len(self.mu)))
        self.P3 = np.zeros((len(self.wav), len(self.mu)))
        self.P4 = np.zeros((len(self.wav), len(self.mu)))

    def normalize_scattering_matrix(self):

        for iw in range(len(self.nu)):

            norm = interp1d_fast_linlog(self.mu, self.P1[iw, :], 0.)

            self.P1[iw, :] /= norm
            self.P2[iw, :] /= norm
            self.P3[iw, :] /= norm
            self.P4[iw, :] /= norm

    def truncate_scattering_matrix(self, mu_max):
        '''
        Remove forward scattering for mu > mu_max
        '''

        self._sort()

        # Loop over wavelengths and reduce scattering cross section
        for iw in range(len(self.wav)):

            # Find fraction remaining
            frac = integrate_linlog_subset(self.mu, self.P1[iw, :],
                                           self.mu[0], mu_max) \
                 / integrate_linlog_subset(self.mu, self.P1[iw, :],
                                           self.mu[0], self.mu[-1])

            logger.info("Removing fraction %g" % frac)

            # Find scattering and absorption opacities
            sigma_nu = self.chi[iw] * self.albedo[iw]
            kappa_nu = self.chi[iw] - sigma

            # Decrease scattering opacity, total opacity, and hence albedo
            sigma_nu *= frac
            self.albedo[iw] = sigma_nu / (sigma_nu + kappa_nu)
            self.chi[iw] = sigma_nu + kappa_nu

        # Interpolate scattering matrix at mu_max
        P1_max = np.zeros((len(self.wav), 1))
        P2_max = np.zeros((len(self.wav), 1))
        P3_max = np.zeros((len(self.wav), 1))
        P4_max = np.zeros((len(self.wav), 1))
        for iw in range(len(self.wav)):
            P1_max[iw, 0] = interp1d_fast_linlog(self.mu, self.P1[iw, :], mu_max)
            P2_max[iw, 0] = interp1d_fast(self.mu, self.P2[iw, :], mu_max)
            P3_max[iw, 0] = interp1d_fast(self.mu, self.P3[iw, :], mu_max)
            P4_max[iw, 0] = interp1d_fast(self.mu, self.P4[iw, :], mu_max)

        # Now truncate scattering matrix elements
        cut = np.searchsorted(self.mu, mu_max)
        self.mu = np.hstack([self.mu[:cut], mu_max])
        self.P1 = np.hstack([self.P1[:, :cut], P1_max])
        self.P2 = np.hstack([self.P2[:, :cut], P2_max])
        self.P3 = np.hstack([self.P3[:, :cut], P3_max])
        self.P4 = np.hstack([self.P4[:, :cut], P4_max])

    def _extrapolate(self, wav1, wav2):

        self._sort()

        nu1 = c / max(wav1, wav2) * 1.e4
        nu2 = c / min(wav1, wav2) * 1.e4

        # need to check both are out of range

        ex_c = extrap1d_log10(self.nu, self.chi)

        self.albedo = np.hstack([self.albedo[0], self.albedo, self.albedo[-1]])
        self.chi = np.hstack([ex_c(nu1), self.chi, ex_c(nu2)])

        self.P1 = np.vstack([self.P1[0, :], self.P1, self.P1[-1, :]])
        self.P2 = np.vstack([self.P2[0, :], self.P2, self.P2[-1, :]])
        self.P3 = np.vstack([self.P3[0, :], self.P3, self.P3[-1, :]])
        self.P4 = np.vstack([self.P4[0, :], self.P4, self.P4[-1, :]])

        self.nu = np.hstack([nu1, self.nu, nu2])
        self.wav = c / self.nu * 1.e4

    def to_table_set(self, table_set):

        # Create optical properties table
        topt = atpy.Table(name='Optical properties')
        topt.add_column('nu', self.nu)
        topt.add_column('albedo', self.albedo)
        topt.add_column('chi', self.chi)

        self.normalize_scattering_matrix()

        topt.add_column('P1', self.P1)
        topt.add_column('P2', self.P2)
        topt.add_column('P3', self.P3)
        topt.add_column('P4', self.P4)

        # Sort by frequency
        topt.sort('nu')

        # Create scattering angles table and add to table set
        tmu = atpy.Table(name='Scattering angles')
        tmu.add_column('mu', self.mu)

        # Add to table set
        table_set.append(topt)
        table_set.append(tmu)

    def from_table_set(self, table_set):

        # Read in the optical properties
        topt = table_set['Optical properties']

        self.nu = topt['nu']
        self.albedo = topt['albedo']
        self.chi = topt['chi']

        self.P1 = topt['P1']
        self.P2 = topt['P2']
        self.P3 = topt['P3']
        self.P4 = topt['P4']

        # Read in the scattering angles
        tmu = table_set['Scattering angles']
        self.mu = tmu['mu']

    def interp_chi_wav(self, wav):
        "Interpolate the opacity to extinction to a given wavelength"
        return interp1d_fast_loglog(self.nu, self.chi, c / (wav * 1.e-4))

    def interp_kappa_wav(self, wav):
        "Interpolate the opacity to absorption to a given wavelength"
        return interp1d_fast_loglog(self.nu, self.kappa, c / (wav * 1.e-4))

    def interp_chi_nu(self, nu):
        "Interpolate the opacity to extinction to a given wavelength"
        return interp1d_fast_loglog(self.nu, self.chi, nu)

    def interp_kappa_nu(self, nu):
        "Interpolate the opacity to absorption to a given wavelength"
        return interp1d_fast_loglog(self.nu, self.kappa, nu)

    def plot(self, figure, subplots):

        import matplotlib.pyplot as plt

        self._sort()

        ax = figure.add_subplot(subplots[0])
        ax.loglog(self.wav, self.chi, color='blue')
        ax2 = ax.twinx()
        ax2.plot(self.wav, self.albedo, color='red')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_ylim(0., 1.)

        self.normalize_scattering_matrix()

        m = plt.cm.gist_heat

        ax = figure.add_subplot(subplots[1])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.abs(self.P1.swapaxes(0, 1))),
                     np.linspace(-2., 2., 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S11', y=0.9, verticalalignment='top', color='white')

        ax = figure.add_subplot(subplots[2])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.abs(self.P2.swapaxes(0, 1))),
                     np.linspace(-2., 2., 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S12', y=0.9, verticalalignment='top', color='white')

        ax = figure.add_subplot(subplots[3])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.abs(self.P3.swapaxes(0, 1))),
                     np.linspace(-2., 2., 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S33', y=0.9, verticalalignment='top', color='white')

        ax = figure.add_subplot(subplots[4])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.abs(self.P4.swapaxes(0, 1))),
                     np.linspace(-2., 2., 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S34', y=0.9, verticalalignment='top', color='white')

        return figure

    def chi_planck_spectrum(self, nu, fnu):
        "Find the Planck mean opacity to extinction for a spectrum"
        if nu.max() > self.nu.max() or nu.min() < self.nu.min():
            raise Exception("Opacity to extinction is not defined at all "
                            "spectrum frequencies")
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        return integrate_loglog(nu, fnu * chi_nu) \
             / integrate_loglog(nu, fnu)

    def kappa_planck_spectrum(self, nu, fnu):
        "Find the Planck mean opacity to absorption for a spectrum"
        if nu.max() > self.nu.max() or nu.min() < self.nu.min():
            raise Exception("Opacity to absorption is not defined at all "
                            "spectrum frequencies")
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        albedo_nu = interp1d_fast_loglog(self.nu, self.albedo, nu)
        kappa_nu = chi_nu * (1. - albedo_nu)
        return integrate_loglog(nu, fnu * kappa_nu) \
             / integrate_loglog(nu, fnu)

    def chi_rosseland_spectrum(self, nu, fnu):
        "Find the Rosseland mean opacity to extinction for a spectrum"
        if nu.max() > self.nu.max() or nu.min() < self.nu.min():
            raise Exception("Opacity to extinction is not defined at all "
                            "spectrum frequencies")
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        return integrate_loglog(nu, fnu) \
             / integrate_loglog(nu, fnu / chi_nu)

    def kappa_rosseland_spectrum(self, nu, fnu):
        "Find the Rosseland mean opacity to absorption for a spectrum"
        if nu.max() > self.nu.max() or nu.min() < self.nu.min():
            raise Exception("Opacity to absorption is not defined at all "
                            "spectrum frequencies")
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        albedo_nu = interp1d_fast_loglog(self.nu, self.albedo, nu)
        kappa_nu = chi_nu * (1. - albedo_nu)
        return integrate_loglog(nu, fnu) \
             / integrate_loglog(nu, fnu / kappa_nu)

    def chi_planck_temperature(self, temperature):
        "Find the Planck mean opacity to extinction for a temperature"

        # Find range of frequencies for Planck function
        planck_nu = planck_nu_range(temperature)
        nu = nu_common(planck_nu, self.nu)

        if planck_nu.min() < self.nu.min():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu >= self.nu.min()]

        if planck_nu.max() > self.nu.max():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu <= self.nu.max()]

        return self.chi_planck_spectrum(nu, B_nu(nu, temperature))

    def kappa_planck_temperature(self, temperature):
        "Find the Rosseland mean opacity to aborption for a temperature"

        # Find range of frequencies for Planck function
        planck_nu = planck_nu_range(temperature)
        nu = nu_common(planck_nu, self.nu)

        if planck_nu.min() < self.nu.min():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu >= self.nu.min()]

        if planck_nu.max() > self.nu.max():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu <= self.nu.max()]

        return self.kappa_planck_spectrum(nu, B_nu(nu, temperature))

    def chi_rosseland_temperature(self, temperature):
        "Find the Rosseland mean opacity to extinction for a temperature"

        # Find range of frequencies for Planck function
        planck_nu = planck_nu_range(temperature)
        nu = nu_common(planck_nu, self.nu)

        if planck_nu.min() < self.nu.min():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu >= self.nu.min()]

        if planck_nu.max() > self.nu.max():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu <= self.nu.max()]

        return self.chi_rosseland_spectrum(nu, B_nu(nu, temperature))

    def kappa_rosseland_temperature(self, temperature):
        "Find the Rosseland mean opacity to absorption for a temperature"

        # Find range of frequencies for Planck function
        planck_nu = planck_nu_range(temperature)
        nu = nu_common(planck_nu, self.nu)

        if planck_nu.min() < self.nu.min():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu >= self.nu.min()]

        if planck_nu.max() > self.nu.max():
            logger.warn("Planck function for requested temperature not completely covered by opacity function")
            nu = nu[nu <= self.nu.max()]

        return self.kappa_rosseland_spectrum(nu, B_nu(nu, temperature))

    def _temperature2specific_energy(self, temperature):

        kappa_planck = self.kappa_planck_temperature(temperature)

        return 4. * sigma * temperature ** 4. * kappa_planck

    def __setattr__(self, attribute, value):
        if attribute in ['nu', 'chi', 'albedo', 'mu'] and value is not None:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError(attribute + " should be a 1-D sequence")
        if attribute in ['nu', 'mu'] and value is not None:
            if not monotonically_increasing(value):
                raise ValueError(attribute + " should be monotonically increasing")
        if attribute == 'nu' and value is not None:
            if value[0] <= 0.:
                raise ValueError('nu should be strictly positive')
        if attribute == 'chi' and value is not None:
            if value[0] < 0.:
                raise ValueError('chi should be positive')
        if attribute == 'albedo' and value is not None:
            if value[0] < 0. or value[-1] > 1.:
                raise ValueError('albedo should be in the range [0:1]')
        if attribute == 'mu' and value is not None:
            if value[0] < -1. or value[-1] > 1.:
                raise ValueError('mu should be in the range [-1:1]')
        if attribute in ['P1', 'P2', 'P3', 'P4'] and value is not None:
            if self.nu is None:
                raise ValueError("nu needs to be set before " + attribute)
            if self.mu is None:
                raise ValueError("mu needs to be set before " + attribute)
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 2:
                raise ValueError(attribute + " should be a 2-D array")
            if value.shape[0] != len(self.nu) or value.shape[1] != len(self.mu):
                raise ValueError(attribute + " has an incorrect shape: %s but expected (%i, %i)" % (value.shape, len(self.nu), len(self.mu)))
        FreezableClass.__setattr__(self, attribute, value)
