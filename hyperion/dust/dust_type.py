import os
import hashlib
import warnings

import atpy
import numpy as np

from hyperion.util.interpolate import interp1d_fast, interp1d_fast_loglog, interp1d_log10, interp1d_fast_linlog
from hyperion.util.integrate import integrate_loglog, integrate_linlog_subset
from hyperion.util.constants import c, sigma
from hyperion.util.functions import B_nu, extrap1d_log10
from hyperion.util.functions import FreezableClass

import matplotlib.pyplot as mpl


def henyey_greenstein(mu, g, p_lin_max):
    P1 = (1.-g*g)/(1.+g*g-2.*g*mu)**1.5
    P2 = - p_lin_max * P1 * (1.-mu*mu)/(1.+mu*mu)
    P3 = P1 * 2. * mu/(1.+mu*mu)
    P4 = 0.
    return P1, P2, P3, P4

class SphericalDust(FreezableClass):

    def __init__(self, *args):

        self.filename = None
        self.md5 = None

        self.n_wav = None
        self.wav = None
        self.nu = None
        self.n_mu = None
        self.mu = None
        self.chi = None
        self.albedo = None
        self.g = None
        self.P1 = None
        self.P2 = None
        self.P3 = None
        self.P4 = None

        self.set_minimum_temperature(0.1)
        self.set_sublimation('no')

        # Set temperature grid for emissivities and mean opacities
        self.n_temp = 1200
        temp_min = 0.1
        temp_max = 100000.
        self.temperatures = np.logspace(np.log10(temp_min), np.log10(temp_max), self.n_temp)
        self.emissivities = False

        self.chi_planck = None
        self.chi_rosseland = None
        self.kappa_planck = None
        self.kappa_rosseland = None

        self.is_lte = True

        self._freeze()

        if len(args) == 0:
            pass
        elif len(args) == 1:
            self.read(args[0])
        else:
            raise Exception("SphericalDust cannot take more than one argument")

    def _temperature2specific_energy_abs(self, temperature):
        return 4. * sigma * temperature ** 4. \
               * self.kappa_planck_temperature(temperature)

    def __getattr__(self, attribute):
        if attribute == 'specific_energy_abs':
            return self._temperature2specific_energy_abs(self.temperatures)
        elif attribute == 'kappa':
            return self.chi * (1. - self.albedo)
        else:
            raise AttributeError("%r object has no attribute %r" %
                                     (type(self).__name__, attribute))

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
            self.g = self.g[::-1]
            self.P1 = self.P1[::-1, :]
            self.P2 = self.P2[::-1, :]
            self.P3 = self.P3[::-1, :]
            self.P4 = self.P4[::-1, :]

    def _extrapolate(self, wav1, wav2):

        self._sort()

        nu1 = c / max(wav1, wav2) * 1.e4
        nu2 = c / min(wav1, wav2) * 1.e4

        # need to check both are out of range

        ex_c = extrap1d_log10(self.nu, self.chi)

        self.albedo = np.hstack([self.albedo[0], self.albedo, self.albedo[-1]])
        self.chi = np.hstack([ex_c(nu1), self.chi, ex_c(nu2)])
        self.g = np.hstack([self.g[0], self.g, self.g[-1]])

        z = np.zeros(self.n_mu)
        self.P1 = np.vstack([self.P1[0, :], self.P1, self.P1[-1, :]])
        self.P2 = np.vstack([self.P2[0, :], self.P2, self.P2[-1, :]])
        self.P3 = np.vstack([self.P3[0, :], self.P3, self.P3[-1, :]])
        self.P4 = np.vstack([self.P4[0, :], self.P4, self.P4[-1, :]])

        self.nu = np.hstack([nu1, self.nu, nu2])
        self.wav = c / self.nu * 1.e4

    def plot(self, filename):

        self._sort()

        mpl.rc('axes', titlesize='x-small')
        mpl.rc('axes', labelsize='x-small')
        mpl.rc('xtick', labelsize='xx-small')
        mpl.rc('ytick', labelsize='xx-small')
        mpl.rc('axes', linewidth=0.5)
        mpl.rc('patch', linewidth=0.5)

        fig = mpl.figure()

        ax1 = fig.add_subplot(3, 2, 1)
        ax1.loglog(self.wav, self.chi, color='blue')
        ax2 = ax1.twinx()
        ax2.plot(self.wav, self.albedo, color='red')
        ax1.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_xlim(self.wav.min(), self.wav.max())
        ax2.set_ylim(0., 1.)

        ax1 = fig.add_subplot(3, 2, 2)

        self.compute_mean_opacities()

        ax1.loglog(self.temperatures, self.chi_planck, color='red')
        ax1.loglog(self.temperatures, self.kappa_planck, color='orange')
        ax1.loglog(self.temperatures, self.chi_rosseland, color='blue')
        ax1.loglog(self.temperatures, self.kappa_rosseland, color='lightblue')
        m = mpl.cm.gist_heat

        self._normalize_scattering_matrix()

        ax1 = fig.add_subplot(3, 2, 3)
        ax1.patch.set_facecolor('black')
        ax1.contourf(self.wav, self.mu, np.log10(np.abs(self.P1.swapaxes(0, 1))), np.linspace(-2., 2., 30), cmap=m)
        ax1.set_xscale('log')
        ax1.set_xlim(self.wav.min(), self.wav.max())
        ax1.set_ylim(-1., 1.)
        ax1.set_title('S11', y=0.9, verticalalignment='top', color='white')

        ax1 = fig.add_subplot(3, 2, 4)
        ax1.patch.set_facecolor('black')
        ax1.contourf(self.wav, self.mu, np.log10(np.abs(self.P2.swapaxes(0, 1))), np.linspace(-2., 2., 30), cmap=m)
        ax1.set_xscale('log')
        ax1.set_xlim(self.wav.min(), self.wav.max())
        ax1.set_ylim(-1., 1.)
        ax1.set_title('S12', y=0.9, verticalalignment='top', color='white')

        ax1 = fig.add_subplot(3, 2, 5)
        ax1.patch.set_facecolor('black')
        ax1.contourf(self.wav, self.mu, np.log10(np.abs(self.P3.swapaxes(0, 1))), np.linspace(-2., 2., 30), cmap=m)
        ax1.set_xscale('log')
        ax1.set_xlim(self.wav.min(), self.wav.max())
        ax1.set_ylim(-1., 1.)
        ax1.set_title('S33', y=0.9, verticalalignment='top', color='white')

        ax1 = fig.add_subplot(3, 2, 6)
        ax1.patch.set_facecolor('black')
        ax1.contourf(self.wav, self.mu, np.log10(np.abs(self.P4.swapaxes(0, 1))), np.linspace(-2., 2., 30), cmap=m)
        ax1.set_xscale('log')
        ax1.set_xlim(self.wav.min(), self.wav.max())
        ax1.set_ylim(-1., 1.)
        ax1.set_title('S34', y=0.9, verticalalignment='top', color='white')

        fig.savefig(filename)

    def _normalize_scattering_matrix(self):

        for iw in range(self.n_wav):

            norm = interp1d_fast_linlog(self.mu, self.P1[iw, :], 0.)

            self.P1[iw, :] /= norm
            self.P2[iw, :] /= norm
            self.P3[iw, :] /= norm
            self.P4[iw, :] /= norm

    def _initialize_scattering_matrix(self):

        self.P1 = np.zeros((self.n_wav, self.n_mu))
        self.P2 = np.zeros((self.n_wav, self.n_mu))
        self.P3 = np.zeros((self.n_wav, self.n_mu))
        self.P4 = np.zeros((self.n_wav, self.n_mu))

    def _truncate_scattering_matrix(self, mu_max):
        '''
        Remove forward scattering for mu > mu_max
        '''

        self._sort()

        # Loop over wavelengths and reduce scattering cross section
        for iw in range(self.n_wav):

            # Find fraction remaining
            frac = integrate_linlog_subset(self.mu, self.P1[iw, :], self.mu[0], mu_max) \
                 / integrate_linlog_subset(self.mu, self.P1[iw, :], self.mu[0], self.mu[-1]) \

            print "Removing fraction %g" % frac

            # Find scattering and absorption opacities
            sigma = self.chi[iw] * self.albedo[iw]
            kappa = self.chi[iw] - sigma

            # Decrease scattering opacity, total opacity, and hence albedo
            sigma *= frac
            self.albedo[iw] = sigma / (sigma + kappa)
            self.chi[iw] = sigma + kappa

        # Interpolate scattering matrix at mu_max
        P1_max = np.zeros((self.n_wav, 1))
        P2_max = np.zeros((self.n_wav, 1))
        P3_max = np.zeros((self.n_wav, 1))
        P4_max = np.zeros((self.n_wav, 1))
        for iw in range(self.n_wav):
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

    def _nu_common(self):

        # Create frequency array for blackbody
        nu_bb_min = self.nu.min()
        nu_bb_max = self.nu.max()
        n_nu_bb = 100
        nu_bb = np.logspace(np.log10(nu_bb_min), np.log10(nu_bb_max), n_nu_bb)
        nu_bb[0] = nu_bb_min
        nu_bb[-1] = nu_bb_max

        # TODO: Ensure that frequency range always contains the BB peak, otherwise mean opacities might be wrong

        # Ensure frequency grid has a minimum spacing
        nu_common = np.array(nu_bb.tolist() + self.nu.tolist())
        nu_common.sort()

        # Remove unique elements (can't just use np.unique because also want to remove very close values)
        keep = (nu_common[1:] - nu_common[:-1])/nu_common[:-1] > 1.e-10
        keep = np.hstack([keep, True])
        nu_common = nu_common[keep]

        return nu_common

    def compute_mean_opacities(self):

        # Find frequencies to use
        nu = self._nu_common()

        # Interpolate opacity to new frequency grid
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        albedo_nu = interp1d_fast_loglog(self.nu, self.albedo, nu)
        kappa_nu = chi_nu * (1. - albedo_nu)

        self.chi_planck = np.zeros(self.n_temp)
        self.kappa_planck = np.zeros(self.n_temp)
        self.chi_rosseland = np.zeros(self.n_temp)
        self.kappa_rosseland = np.zeros(self.n_temp)

        # Loop through the temperatures and compute mean opacities
        for it, T in enumerate(self.temperatures):

            # Compute planck mean opacity
            self.chi_planck[it] = integrate_loglog(nu, B_nu(nu, T) * chi_nu) \
                               / integrate_loglog(nu, B_nu(nu, T))

            # Compute planck mean absoptive opacity
            self.kappa_planck[it] = integrate_loglog(nu, B_nu(nu, T) * kappa_nu) \
                                 / integrate_loglog(nu, B_nu(nu, T))

            # Compute Rosseland mean opacity
            self.chi_rosseland[it] = integrate_loglog(nu, B_nu(nu, T)) \
                                   / integrate_loglog(nu, B_nu(nu, T) / chi_nu)

            # Compute Rosseland mean opacity
            self.kappa_rosseland[it] = integrate_loglog(nu, B_nu(nu, T)) \
                                     / integrate_loglog(nu, B_nu(nu, T) / kappa_nu)

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

    def chi_planck_temperature(self, temperature):
        "Find the Planck mean opacity to extinction for a blackbody spectrum"
        return interp1d_fast_loglog(self.temperatures, self.chi_planck, temperature)

    def kappa_planck_temperature(self, temperature):
        "Find the Planck mean opacity to extinction for a blackbody spectrum"
        return interp1d_fast_loglog(self.temperatures, self.kappa_planck, temperature)

    def chi_rosseland_temperature(self, temperature):
        "Find the Rosseland mean opacity to extinction for a blackbody spectrum"
        return interp1d_fast_loglog(self.temperatures, self.chi_rosseland, temperature)

    def kappa_rosseland_temperature(self, temperature):
        "Find the Rosseland mean opacity to extinction for a blackbody spectrum"
        return interp1d_fast_loglog(self.temperatures, self.kappa_rosseland, temperature)

    def chi_planck_spectrum(self, nu, fnu):
        "Find the Planck mean opacity to extinction for a spectrum"
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        return integrate_loglog(nu, fnu * chi_nu) \
             / integrate_loglog(nu, fnu)

    def kappa_planck_spectrum(self, nu, fnu):
        "Find the Planck mean opacity to absorption for a spectrum"
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        albedo_nu = interp1d_fast_loglog(self.nu, self.albedo, nu)
        kappa_nu = chi_nu * (1. - albedo_nu)
        return integrate_loglog(nu, fnu * kappa_nu) \
             / integrate_loglog(nu, fnu)

    def chi_rosseland_spectrum(self, nu, fnu):
        "Find the Rosseland mean opacity to extinction for a spectrum"
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        return integrate_loglog(nu, fnu) \
             / integrate_loglog(nu, fnu / chi_nu)

    def kappa_rosseland_spectrum(self, nu, fnu):
        "Find the Rosseland mean opacity to absorption for a spectrum"
        chi_nu = interp1d_fast_loglog(self.nu, self.chi, nu)
        albedo_nu = interp1d_fast_loglog(self.nu, self.albedo, nu)
        kappa_nu = chi_nu * (1. - albedo_nu)
        return integrate_loglog(nu, fnu) \
             / integrate_loglog(nu, fnu / kappa_nu)

    def set_sublimation(self, mode, temperature=None):
        '''
        Set the dust sublimation parameters.

        Parameters
        ----------
        mode : str
            The dust sublimation mode, which can be:
                * 'no'   - no sublimation
                * 'fast' - remove all dust in cells exceeding the
                           sublimation temperature/energy
                * 'slow' - reduce the dust in cells exceeding the
                           sublimation temperature/energy
                * 'cap'  - any temperature/energy exceeding the sublimation
                           value is reset to the sublimation value.

        temperature : float, optional
            The dust sublimation temperature, in K
        '''

        if mode not in ['no', 'fast', 'slow', 'cap']:
            raise Exception("mode should be one of no/fast/slow/cap")

        if mode != 'no' and temperature is None:
            raise Exception("Need to specify a sublimation temperature")

        self.sublimation_mode = mode
        self.sublimation_temperature = temperature

    def set_minimum_temperature(self, temperature):
        '''
        Set the minimum dust temperature

        Dust which has a temperature that falls below this value will be
        reset to the minimum at the end of each iteration.

        Parameters
        ----------
        temperature : float
            The minimum temperature in K
        '''
        self.minimum_temperature = temperature

    def write(self, filename, compression=True):
        '''
        Write out to a standard dust file, including calculations of the mean
        opacities and optionally thermal emissivities.
        '''

        self._sort()

        # Create dust table set
        ts = atpy.TableSet()

        # Add standard keywords to header
        ts.add_keyword('version', 1)
        ts.add_keyword('type', 1)
        if self.md5:
            ts.add_keyword('asciimd5', self.md5)

        # Create optical properties table
        topt = atpy.Table(name='Optical properties')
        topt.add_column('nu', self.nu)
        topt.add_column('albedo', self.albedo)
        topt.add_column('chi', self.chi)

        self._normalize_scattering_matrix()

        topt.add_column('P1', self.P1)
        topt.add_column('P2', self.P2)
        topt.add_column('P3', self.P3)
        topt.add_column('P4', self.P4)

        # Sort by frequency
        topt.sort('nu')

        # Add to table set
        ts.append(topt)

        # Create scattering angles table and add to table set
        tmu = atpy.Table(name='Scattering angles')
        tmu.add_column('mu', self.mu)
        ts.append(tmu)

        # Create table to contain the planck and rosseland mean opacities

        self.compute_mean_opacities()

        tmean = atpy.Table(name='Mean opacities')

        tmean.add_column('specific_energy_abs', self.specific_energy_abs)
        tmean.add_column('temperature', self.temperatures)

        tmean.add_column('chi_planck', self.chi_planck)
        tmean.add_column('kappa_planck', self.kappa_planck)
        tmean.add_column('chi_rosseland', self.chi_rosseland)
        tmean.add_column('kappa_rosseland', self.kappa_rosseland)

        # Add to table set
        ts.append(tmean)

        if self.emissivities:

            # Assume emissivities are not LTE
            self.is_lte = False

            # Read in existing emissivities
            te = atpy.Table(self.emissivities, table='Emissivities')

            # Set frequency scale
            emiss_nu = te.nu

            # Set emissivities
            jnu = te.jnu

            # Read in emissivity variable
            tev = atpy.Table(self.emissivities, table='Emissivity variable')

            # Set emissivity variable
            emissvar = tev.names[0]

            # Only specific energy is considered a valid emissivity variable
            if tev.names[0] != 'specific_energy_abs':
                raise Exception("Unknown emissivity variable: %s" % emissvar)

            # Find specific energy emissivity variable
            specific_energy_abs = tev.specific_energy_abs

        else: # Compute LTE emissivities

            # Emissivities are LTE
            self.is_lte = True

            # Set frequency scale
            emiss_nu = self._nu_common()

            # Compute opacity to absorption
            kappa_nu = interp1d_fast_loglog(self.nu, self.kappa, emiss_nu)

            # Compute LTE emissivities
            jnu = np.zeros((len(emiss_nu), self.n_temp))
            for it, T in enumerate(self.temperatures):
                jnu[:, it] = kappa_nu * B_nu(emiss_nu, T)

            # Set emissivity variable to specific energy
            emissvar = 'specific_energy_abs'

            # Compute specific energy emissivity variable
            specific_energy_abs = self.specific_energy_abs

        # Add header keyword to specify emissivity variable
        if emissvar == 'specific_energy_abs':
            ts.add_keyword('emissvar', 'E')
        else:
            raise Exception("Unknown emissivity variable: %s" % emissvar)

        # Add header keyword to specify whether dust is LTE
        ts.add_keyword('lte', 'yes' if self.is_lte else 'no')

        # Add emissivities to table set
        temiss = atpy.Table(name='Emissivities')
        temiss.add_column('nu', emiss_nu)
        temiss.add_column('jnu', jnu)
        ts.append(temiss)

        # Add emissivity variable to table set
        temissvar = atpy.Table(name='Emissivity variable')
        if emissvar == 'specific_energy_abs':
            temissvar.add_column('specific_energy_abs', specific_energy_abs)
        else:
            raise Exception("Unknown emissivity variable: %s" % emissvar)
        ts.append(temissvar)

        # Dust sublimation parameters
        ts.add_keyword('sublimation_mode', self.sublimation_mode)
        if self.sublimation_mode != 'no':
            ts.add_keyword('sublimation_specific_energy', self._temperature2specific_energy_abs(self.sublimation_temperature))

        # Minimum temperature parameter
        ts.add_keyword('minimum_specific_energy', self._temperature2specific_energy_abs(self.minimum_temperature))

        # Output dust file
        ts.write(filename, overwrite=True, compression=compression, type='hdf5')

        self.filename = filename

    def read(self, filename):
        '''
        Read in from a standard dust file
        '''

        # Check file exists
        if not os.path.exists(filename):
            raise Exception("File not found: %s" % filename)

        self.filename = filename

        # Read in dust table set
        ts = atpy.TableSet(filename, verbose=False)

        # Check version and type
        if ts.keywords['version'] != 1:
            raise Exception("Version should be 1")
        if ts.keywords['type'] != 1:
            raise Exception("Type should be 1")
        if 'asciimd5' in ts.keywords:
            self.md5 = ts.keywords['asciimd5']
        else:
            self.md5 = None

        # Read in the optical properties
        topt = ts['Optical properties']
        self.n_wav = len(topt)
        self.nu = topt['nu']
        self.albedo = topt['albedo']
        self.chi = topt['chi']

        self.P1 = topt['P1']
        self.P2 = topt['P2']
        self.P3 = topt['P3']
        self.P4 = topt['P4']

        # Read in the scattering angles
        tmu = ts['Scattering angles']
        self.mu = tmu['mu']

        # Read in the planck and rosseland mean opacities
        tmean = ts['Mean opacities']

        self.temperatures = tmean['temperature']

        self.chi_planck = tmean['chi_planck']
        self.kappa_planck = tmean['kappa_planck']
        self.chi_rosseland = tmean['chi_rosseland']
        self.kappa_rosseland = tmean['kappa_rosseland']

        # Find the emissivity variable type
        if ts.keywords['emissvar'] == 'E':
            emissvar = 'specific_energy_abs'
        else:
            raise Exception("Unknown emissivity variable: %s" % ts.keywords['emissvar'])

        # Read emissivities
        temiss = ts['Emissivities']
        emiss_nu = temiss['nu']
        jnu = temiss['jnu']

        # Read in emissivity variable
        temissvar = ts['Emissivity variable']
        if emissvar == 'specific_energy_abs':
            specific_energy_abs = temissvar['specific_energy_abs']
        else:
            raise Exception("Unknown emissivity variable: %s" % emissvar)


class IsotropicSphericalDust(SphericalDust):

    def __init__(self, wav, chi, albedo):

        SphericalDust.__init__(self)

        # Set cos(theta) grid for computing the scattering matrix elements
        self.n_mu = 2
        self.mu = np.linspace(-1., 1., self.n_mu)

        self.n_wav = len(chi)
        self.wav = wav
        self.nu = c / self.wav * 1.e4
        self.albedo = albedo
        self.chi = chi
        self.g = np.ones(self.n_wav)

        # Compute scattering matrix elements
        self._initialize_scattering_matrix()

        self.P1[:, :] = 1.
        self.P2[:, :] = 0.
        self.P3[:, :] = 1.
        self.P4[:, :] = 0.


class SimpleSphericalDust(SphericalDust):

    def __init__(self, filename):

        SphericalDust.__init__(self)

        # Set cos(theta) grid for computing the scattering matrix elements
        self.n_mu = 100
        self.mu = np.linspace(-1., 1., self.n_mu)

        # Read in dust file
        dustfile = np.loadtxt(filename, dtype=[('wav', float), ('c_ext', float), \
                              ('c_sca', float), ('chi', float), ('g', float), \
                              ('p_lin_max', float)], usecols=[0, 1, 2, 3, 4, 5])

        self.n_wav = len(dustfile)
        self.wav = dustfile['wav']
        self.nu = c / self.wav * 1.e4
        self.albedo = dustfile['c_sca'] / dustfile['c_ext']
        self.chi = dustfile['chi']
        self.g = dustfile['g']

        self.md5 = hashlib.md5(open(filename, 'rb').read()).hexdigest()

        # Compute scattering matrix elements
        self._initialize_scattering_matrix()

        for i in range(0, self.n_mu):
            self.P1[:, i], self.P2[:, i], self.P3[:,i], self.P4[:,i] = henyey_greenstein(self.mu[i], self.g, dustfile['p_lin_max'])

class CoatsphSingle(SphericalDust):

    def __init__(self, directory, size, density):
        '''
        Initialize single-component dust.

        Required Arguments:

            *directory*: [ string ]
                Directory containing all the files describing the dust

            *size*: [ float ]
                Grain size, in cm

            *density*: [ float ]
                Dust grain density, in g/cm^3
        '''

        SphericalDust.__init__(self)

        f = open('%s/coatsph_forw.dat' % directory, 'rb')
        version = f.readline()
        n_components = int(f.readline().strip().split()[5])

        # Read in main dust file

        dustfile = np.loadtxt(f, skiprows=3,
                    dtype=[('x', float), ('radius', float), ('wav', float),
                    ('q_ext', float), ('q_sca', float), ('q_back', float),
                    ('g', float)])

        self.n_wav = len(dustfile)

        self.wav = np.zeros(self.n_wav)
        self.chi = np.zeros(self.n_wav)
        self.albedo = np.zeros(self.n_wav)

        self.wav = dustfile['wav']
        self.nu = c / self.wav * 1.e4
        self.albedo = dustfile['q_sca'] / dustfile['q_ext']
        self.chi = 0.75 * dustfile['q_ext'] / size / density
        self.g = dustfile['g']

        # Read in scattering matrix elements

        for i in range(self.n_wav):

            filename = '%s/coatsph_scat_%04i_0001.dat' % (directory, i+1)

            phasefile = np.loadtxt(filename, skiprows=9,
                        dtype=[('theta', float), ('s11', float), ('polariz',
                        float), ('s12', float), ('s33', float), ('s34',
                        float)])

            if i==0:
                self.n_mu = len(phasefile)
                self.mu = np.cos(np.radians(phasefile['theta']))
                self._initialize_scattering_matrix()

            self.P1[i, :] = phasefile['s11']
            self.P2[i, :] = phasefile['s12']
            self.P3[i, :] = phasefile['s33']
            self.P4[i, :] = phasefile['s34']


class CoatsphMultiple(SphericalDust):

    def __init__(self, directory):
        '''
        Initialize multi-component dust.

        Required Arguments:

            *directory*: [ string ]
                Directory containing all the files describing the dust
        '''

        SphericalDust.__init__(self)

        f = open('%s/coatsph_forw.dat' % directory, 'rb')
        version = f.readline()
        n_components = int(f.readline().strip().split()[5])

        # Read in main dust file

        dustfile = np.loadtxt(f, skiprows=7,
                    dtype=[('wav', float), ('c_ext', float), ('c_sca', float),
                    ('chi', float), ('g', float), ('pmax', float),
                    ('thetmax', float)])

        self.n_wav = len(dustfile)
        self.wav = dustfile['wav']
        self.nu = c / self.wav * 1.e4
        self.albedo = dustfile['c_sca'] / dustfile['c_ext']
        self.chi = dustfile['chi']
        self.g = dustfile['g']

        # Read in scattering matrix elements

        for i in range(self.n_wav):

            filename = '%s/coatsph_scat.%04i.dat' % (directory, i+1)

            phasefile = np.loadtxt(filename, skiprows=7,
                        dtype=[('theta', float), ('s11', float), ('polariz',
                        float), ('s12', float), ('s33', float), ('s34',
                        float)])

            if i==0:
                self.n_mu = len(phasefile)
                self.mu = np.cos(np.radians(phasefile['theta']))
                self._initialize_scattering_matrix()

            self.P1[i, :] = phasefile['s11']
            self.P2[i, :] = phasefile['s12']
            self.P3[i, :] = phasefile['s33']
            self.P4[i, :] = phasefile['s34']


class MieXDust(SphericalDust):

    def __init__(self, model):

        SphericalDust.__init__(self)

        self.wav = np.loadtxt('%s.alb' % model, usecols=[0])
        self.albedo = np.loadtxt('%s.alb' % model, usecols=[1])
        self.g = np.loadtxt('%s.g' % model, usecols=[1])
        kappa = np.loadtxt('%s.k_abs' % model, usecols=[1])
        self.chi = kappa / (1 - self.albedo)

        # Check for NaN values
        for quantity in ['chi', 'albedo']:

            values = self.__dict__[quantity]

            if np.any(np.isnan(values)):
                warnings.warn("NaN values found inside MieX %s file - interpolating" % quantity)
                invalid = np.isnan(values)
                f = interp1d_log10(self.wav[~invalid], values[~invalid])
                values[invalid] = f(self.wav[invalid])
                if np.any(np.isnan(values)):
                    raise Exception("Did not manage to fix NaN values in MieX %s" % quantity)

        self.nu = c / self.wav * 1.e4

        self.n_wav = len(self.wav)
        self.n_mu = (len(open('%s.f11' % model).readlines()) / self.n_wav) - 1

        self.mu = np.zeros(self.n_mu)
        self._initialize_scattering_matrix()

        # Read mu
        f11 = open('%s.f11' % model)
        f11.readline()
        f11.readline()
        for i in range(self.n_mu):
            self.mu[i] = np.cos(np.radians(float(f11.readline().split()[0])))
        f11.close()

        # Read in matrix elements
        f11 = open('%s.f11' % model)
        f12 = open('%s.f12' % model)
        f33 = open('%s.f33' % model)
        f34 = open('%s.f34' % model)

        f11.readline()
        f12.readline()
        f33.readline()
        f34.readline()

        for j in range(self.n_wav):

            if float(f11.readline()) != self.wav[j]:
                raise Exception("Incorrect wavelength in f11")
            if float(f12.readline()) != self.wav[j]:
                raise Exception("Incorrect wavelength in f12")
            if float(f33.readline()) != self.wav[j]:
                raise Exception("Incorrect wavelength in f33")
            if float(f34.readline()) != self.wav[j]:
                raise Exception("Incorrect wavelength in f34")

            for i in range(self.n_mu):

                self.P1[j, i] = float(f11.readline().split()[1])
                self.P2[j, i] = float(f12.readline().split()[1])
                self.P3[j, i] = float(f33.readline().split()[1])
                self.P4[j, i] = float(f34.readline().split()[1])

        for i in range(self.n_mu):

            for quantity in ['P1', 'P2', 'P3', 'P4']:

                values = self.__dict__[quantity]

                if np.any(np.isnan(values[:, i])):
                    warnings.warn("NaN values found inside MieX %s file - interpolating" % quantity)
                    invalid = np.isnan(values[:, i])
                    f = interp1d_log10(self.wav[~invalid], values[:, i][~invalid])
                    values[:, i][invalid] = f(self.wav[invalid])
                    if np.any(np.isnan(values[:, i])):
                        raise Exception("Did not manage to fix NaN values in MieX %s" % quantity)


class BHDust(SphericalDust):

    def __init__(self, model):

        SphericalDust.__init__(self)

        self.wav = np.loadtxt('%s.wav' % model)
        self.mu = np.loadtxt('%s.mu' % model)
        self.albedo = np.loadtxt('%s.alb' % model)
        self.chi = np.loadtxt('%s.chi' % model)
        self.g = np.loadtxt('%s.g' % model)

        self.nu = c / self.wav * 1.e4

        self.n_wav = len(self.wav)
        self.n_mu = len(self.mu)

        self._initialize_scattering_matrix()

        self.P1 = np.loadtxt('%s.f11' % model)
        self.P2 = np.loadtxt('%s.f12' % model)
        self.P3 = np.loadtxt('%s.f33' % model)
        self.P4 = np.loadtxt('%s.f34' % model)
