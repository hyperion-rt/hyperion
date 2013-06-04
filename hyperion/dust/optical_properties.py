from __future__ import print_function, division

import hashlib

import numpy as np
from astropy.table import Table, Column

from ..util.integrate import integrate_loglog, integrate_linlog_subset
from ..util.interpolate import interp1d_fast, interp1d_fast_loglog, \
                                      interp1d_fast_linlog
from ..util.functions import extrap1d_log10, B_nu, FreezableClass, \
                                    nu_common, planck_nu_range, is_numpy_array, monotonically_increasing
from ..util.constants import c, sigma
from astropy import log as logger


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

            self.nu = self.nu[::-1]
            self.albedo = self.albedo[::-1]
            self.chi = self.chi[::-1]
            self.P1 = self.P1[::-1, :]
            self.P2 = self.P2[::-1, :]
            self.P3 = self.P3[::-1, :]
            self.P4 = self.P4[::-1, :]

    def initialize_scattering_matrix(self):

        self.P1 = np.zeros((len(self.nu), len(self.mu)))
        self.P2 = np.zeros((len(self.nu), len(self.mu)))
        self.P3 = np.zeros((len(self.nu), len(self.mu)))
        self.P4 = np.zeros((len(self.nu), len(self.mu)))

    def normalize_scattering_matrix(self):

        for inu in range(len(self.nu)):

            norm = interp1d_fast_linlog(self.mu, self.P1[inu, :], 0.)

            self.P1[inu, :] /= norm
            self.P2[inu, :] /= norm
            self.P3[inu, :] /= norm
            self.P4[inu, :] /= norm

    def truncate_scattering_matrix(self, mu_max):
        '''
        Remove forward scattering for mu > mu_max
        '''

        self._sort()

        # Loop over wavelengths and reduce scattering cross section
        for inu in range(len(self.nu)):

            # Find fraction remaining
            frac = integrate_linlog_subset(self.mu, self.P1[inu, :],
                                           self.mu[0], mu_max) \
                 / integrate_linlog_subset(self.mu, self.P1[inu, :],
                                           self.mu[0], self.mu[-1])

            logger.info("Removing fraction %g" % frac)

            # Find scattering and absorption opacities
            sigma_nu = self.chi[inu] * self.albedo[inu]
            kappa_nu = self.chi[inu] - sigma

            # Decrease scattering opacity, total opacity, and hence albedo
            sigma_nu *= frac
            self.albedo[inu] = sigma_nu / (sigma_nu + kappa_nu)
            self.chi[inu] = sigma_nu + kappa_nu

        # Interpolate scattering matrix at mu_max
        P1_max = np.zeros((len(self.nu), 1))
        P2_max = np.zeros((len(self.nu), 1))
        P3_max = np.zeros((len(self.nu), 1))
        P4_max = np.zeros((len(self.nu), 1))
        for inu in range(len(self.nu)):
            P1_max[inu, 0] = interp1d_fast_linlog(self.mu, self.P1[inu, :], mu_max)
            P2_max[inu, 0] = interp1d_fast(self.mu, self.P2[inu, :], mu_max)
            P3_max[inu, 0] = interp1d_fast(self.mu, self.P3[inu, :], mu_max)
            P4_max[inu, 0] = interp1d_fast(self.mu, self.P4[inu, :], mu_max)

        # Now truncate scattering matrix elements
        cut = np.searchsorted(self.mu, mu_max)
        self.mu = np.hstack([self.mu[:cut], mu_max])
        self.P1 = np.hstack([self.P1[:, :cut], P1_max])
        self.P2 = np.hstack([self.P2[:, :cut], P2_max])
        self.P3 = np.hstack([self.P3[:, :cut], P3_max])
        self.P4 = np.hstack([self.P4[:, :cut], P4_max])

    def extrapolate_wav(self, wav1, wav2):
        '''
        Extrapolate the optical properties to a larger frequency range.

        Parameters
        ----------
        wav1, wav2 : float
            The range of wavelengths (in microns) to extrapolate the optical
            properties to.

        Notes
        -----
        The extrapolation is done in the following way:

            * The opacity to extinction (``chi``) is extrapolated by fitting a
              power-law to the opacities at the two highest frequencies and
              following that power law, and similarly at the lowest
              frequencies. This ensures that the slope of the opacity remains
              constant.

            * The albedo is extrapolated by assuming that the albedo is
              constant outside the original range, and is set to the same
              value as the values for the lowest and highest frequencies.

            * The scattering matrix is extrapolated similarly to the albedo,
              by simply extending the values for the lowest and highest
              frequencies to the new frequency range.
        '''

        nu1 = c / max(wav1, wav2) * 1.e4
        nu2 = c / min(wav1, wav2) * 1.e4

        return self.extrapolate_nu(nu1, nu2)

    def extrapolate_nu(self, nu1, nu2):
        '''
        Extrapolate the optical properties to a larger frequency range.

        Parameters
        ----------
        nu1, nu2 : float
            The range of frequencies to extrapolate the optical properties to.

        Notes
        -----
        The extrapolation is done in the following way:

            * The opacity to extinction (``chi``) is extrapolated by fitting a
              power-law to the opacities at the two highest frequencies and
              following that power law, and similarly at the lowest
              frequencies. This ensures that the slope of the opacity remains
              constant.

            * The albedo is extrapolated by assuming that the albedo is
              constant outside the original range, and is set to the same
              value as the values for the lowest and highest frequencies.

            * The scattering matrix is extrapolated similarly to the albedo,
              by simply extending the values for the lowest and highest
              frequencies to the new frequency range.
        '''

        self._sort()

        if nu1 >= self.nu[0]:

            logger.info("Lower frequency is inside existing range, no extrapolation will be done at the lowest frequencies")

        else:

            ex_c = extrap1d_log10(self.nu, self.chi)

            self.albedo = np.hstack([self.albedo[0], self.albedo])
            self.chi = np.hstack([ex_c(nu1), self.chi])
            self.nu = np.hstack([nu1, self.nu])

            self.P1 = np.vstack([self.P1[0, :], self.P1])
            self.P2 = np.vstack([self.P2[0, :], self.P2])
            self.P3 = np.vstack([self.P3[0, :], self.P3])
            self.P4 = np.vstack([self.P4[0, :], self.P4])

        if nu2 <= self.nu[-1]:

            logger.info("Upper frequency is inside existing range, no extrapolation will be done at the highest frequencies")

        else:

            ex_c = extrap1d_log10(self.nu, self.chi)

            self.albedo = np.hstack([self.albedo, self.albedo[-1]])
            self.chi = np.hstack([self.chi, ex_c(nu2)])
            self.nu = np.hstack([self.nu, nu2])

            self.P1 = np.vstack([self.P1, self.P1[-1, :]])
            self.P2 = np.vstack([self.P2, self.P2[-1, :]])
            self.P3 = np.vstack([self.P3, self.P3[-1, :]])
            self.P4 = np.vstack([self.P4, self.P4[-1, :]])

    def to_hdf5_group(self, group):

        if not self.all_set():
            raise Exception("Not all attributes of the optical properties are set")

        # Create optical properties table
        topt = Table()
        topt.add_column(Column(data=self.nu, name='nu'))
        topt.add_column(Column(data=self.albedo, name='albedo'))
        topt.add_column(Column(data=self.chi, name='chi'))

        self.normalize_scattering_matrix()

        topt.add_column(Column(data=self.P1, name='P1'))
        topt.add_column(Column(data=self.P2, name='P2'))
        topt.add_column(Column(data=self.P3, name='P3'))
        topt.add_column(Column(data=self.P4, name='P4'))

        # Sort by frequency
        topt.sort('nu')

        # Create scattering angles table and add to table set
        tmu = Table()
        tmu.add_column(Column(data=self.mu, name='mu'))

        # Add to group
        topt.write(group, path='optical_properties')
        tmu.write(group, path='scattering_angles')

    def from_hdf5_group(self, group):

        # Read in the scattering angles
        tmu = group['scattering_angles']
        self.mu = tmu['mu']

        # Read in the optical properties
        topt = group['optical_properties']

        self.nu = topt['nu']
        self.albedo = topt['albedo']
        self.chi = topt['chi']

        self.P1 = topt['P1']
        self.P2 = topt['P2']
        self.P3 = topt['P3']
        self.P4 = topt['P4']

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

    def all_set(self):
        return self.get_missing_attributes() == []

    def get_missing_attributes(self):
        missing = []
        for attribute in ['nu', 'chi', 'albedo', 'mu', 'P1', 'P2', 'P3', 'P4']:
            if getattr(self, attribute) is None:
                missing.append(attribute)
        return missing

    def ensure_all_set(self):
        if not self.all_set():
            missing = self.get_missing_attributes()
            raise Exception("The following attributes of the optical properties have not been set: {0:s}".format(', '.join(missing)))

    def plot(self, figure, subplots):

        self.ensure_all_set()

        import matplotlib.pyplot as plt

        self._sort()

        ax = figure.add_subplot(subplots[0])
        ax.loglog(self.wav, self.chi, color='blue')
        ax2 = ax.twinx()
        ax2.plot(self.wav, self.albedo, color='red')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_ylim(0., 1.)
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('Opacity to extinction (cm^2/g)', color='blue')
        ax2.set_ylabel('Albedo', color='red', rotation=-90)
        self.normalize_scattering_matrix()

        m = plt.cm.gist_heat
        vmin, vmax = -2., 2.

        ax = figure.add_subplot(subplots[1])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.clip(np.abs(self.P1.swapaxes(0, 1)), 10. ** vmin, 10. ** vmax)),
                     np.linspace(vmin, vmax, 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S11', y=0.9, verticalalignment='top', color='white')
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('mu')

        ax = figure.add_subplot(subplots[2])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.clip(np.abs(self.P2.swapaxes(0, 1)), 10. ** vmin, 10. ** vmax)),
                     np.linspace(vmin, vmax, 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S12', y=0.9, verticalalignment='top', color='white')
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('mu')

        ax = figure.add_subplot(subplots[3])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.clip(np.abs(self.P3.swapaxes(0, 1)), 10. ** vmin, 10. ** vmax)),
                     np.linspace(vmin, vmax, 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S33', y=0.9, verticalalignment='top', color='white')
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('mu')

        ax = figure.add_subplot(subplots[4])
        ax.patch.set_facecolor('black')
        ax.contourf(self.wav, self.mu,
                     np.log10(np.clip(np.abs(self.P4.swapaxes(0, 1)), 10. ** vmin, 10. ** vmax)),
                     np.linspace(vmin, vmax, 30), cmap=m)
        ax.set_xscale('log')
        ax.set_xlim(self.wav.min(), self.wav.max())
        ax.set_ylim(-1., 1.)
        ax.set_title('S34', y=0.9, verticalalignment='top', color='white')
        ax.set_xlabel('Wavelength (microns)')
        ax.set_ylabel('mu')

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

    def hash(self):
        h = hashlib.md5()
        h.update(self.nu.tostring())
        h.update(self.chi.tostring())
        h.update(self.albedo.tostring())
        h.update(self.mu.tostring())
        h.update(self.P1.tostring())
        h.update(self.P2.tostring())
        h.update(self.P3.tostring())
        h.update(self.P4.tostring())
        return h.hexdigest()

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
