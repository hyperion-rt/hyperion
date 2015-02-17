##############################################################################
#
#  File format versions:
#
#  1: initial version
#
#  2: now contains reciprocal planck opacity and rosseland opacity
#     previous rosseland opacity was actually reciprocal planck opacity
#
##############################################################################

from __future__ import print_function, division

import os
import hashlib

import h5py
import numpy as np
from astropy import log as logger
from astropy.extern import six

from ..version import __version__

from ..util.constants import c
from ..util.functions import FreezableClass
from ..util.interpolate import interp1d_fast_loglog
from ..util.integrate import integrate_loglog

from .optical_properties import OpticalProperties
from .emissivities import Emissivities
from .mean_opacities import MeanOpacities


def henyey_greenstein(mu, g, p_lin_max):
    P1 = (1. - g * g) / (1. + g * g - 2. * g * mu) ** 1.5
    P2 = - p_lin_max * P1 * (1. - mu * mu) / (1. + mu * mu)
    P3 = P1 * 2. * mu / (1. + mu * mu)
    P4 = 0.
    return P1, P2, P3, P4


class SphericalDust(FreezableClass):
    r"""
    This class should be used in cases where fully arbitrary dust properties
    need to be specified, within the framework of randomly oriented grains,
    which means that the scattering phase function has the general form:

    .. math:: R(\theta) = \left[\begin{array}{cccc}P_1 & P_2 & 0 & 0 \\P_2 & P_1 & 0 & 0 \\0 & 0 & P_3 & -P_4 \\0 & 0 & P_4 & P_3\end{array}\right]

    This class is initialized with::

        d = SphericalDust()

    and the properties should then be set manually. See
    `here <http://docs.hyperion-rt.org/en/stable/setup/setup_dust.html#fully-customized-4-element-dust>`_
    for a description of the available properties and how to set them.
    """

    def __init__(self, *args):

        self._file = None
        self.md5 = None

        self.optical_properties = OpticalProperties()
        self.mean_opacities = MeanOpacities()
        self.emissivities = Emissivities()

        self.set_sublimation_specific_energy('no', 0.)

        self._freeze()

        if len(args) == 0:
            pass
        elif len(args) == 1:
            self.read(args[0])
        else:
            raise Exception("SphericalDust cannot take more than one argument")

    def hash(self):

        h = hashlib.md5()

        if self.optical_properties is None or not self.optical_properties.all_set():
            h.update('none'.encode('utf-8'))
        else:
            h.update(self.optical_properties.hash().encode('utf-8'))

        if self.emissivities is None or not self.emissivities.all_set():
            h.update('none'.encode('utf-8'))
        else:
            h.update(self.emissivities.hash().encode('utf-8'))

        if self.mean_opacities is None or not self.mean_opacities.all_set():
            h.update('none'.encode('utf-8'))
        else:
            h.update(self.mean_opacities.hash().encode('utf-8'))

        import struct
        h.update(self.sublimation_mode.encode('utf-8'))
        h.update(struct.pack('>d', self.sublimation_energy))

        return h.hexdigest()

    def set_lte_emissivities(self, n_temp=1200, temp_min=0.1, temp_max=100000.):
        '''
        Calculate the emissivities assuming LTE

        Parameters
        ----------
        n_temp : int, optional
            The number of temperatures to calculate the emissivities for
        temp_min : float, optional
            The minimum temperature to calculate the emissivities for
        temp_max : float, optional
            The maximum temperature to calculate the emissivities for
        '''
        self.mean_opacities.compute(self.optical_properties, n_temp=n_temp,
                                    temp_min=temp_min, temp_max=temp_max)
        self.emissivities.set_lte(self.optical_properties, self.mean_opacities)

    def plot(self, filename):

        # Check that the optical properties have been set
        self.optical_properties.ensure_all_set()

        import matplotlib.pyplot as plt

        # Save original rc parameters
        rc_orig = plt.rcParams

        # Reset to defaults
        plt.rcdefaults()
        plt.rc('legend', fontsize=7)
        plt.rc('axes', titlesize='x-small')
        plt.rc('axes', labelsize='x-small')
        plt.rc('xtick', labelsize='xx-small')
        plt.rc('ytick', labelsize='xx-small')
        plt.rc('axes', linewidth=0.5)
        plt.rc('patch', linewidth=0.5)

        # Compute mean opacities if not already existent
        self._compute_mean_opacities()

        # Check that emissivities are set (before computing mean opacities)
        if not self.emissivities.all_set():
            logger.info("Computing emissivities assuming LTE")
            self.emissivities.set_lte(self.optical_properties, self.mean_opacities)

        # Initialize figure
        fig = plt.figure(figsize=(10, 12))

        # Plot optical properties
        fig = self.optical_properties.plot(fig, [421, 423, 424, 425, 426])

        # Plot mean opacities
        fig = self.mean_opacities.plot(fig, 428)

        # Plot emissivities
        fig = self.emissivities.plot(fig, 427)

        # Adjust spacing between subplots
        fig.subplots_adjust(left=0.08, right=0.92, wspace=0.22, hspace=0.30)

        # Save figure
        fig.savefig(filename, bbox_inches='tight')

        # Close figure to save RAM
        plt.close(fig)

        # Restore rc parameters
        plt.rc(rc_orig)

    def set_sublimation_temperature(self, mode, temperature=0.):
        '''
        Set the dust sublimation mode and temperature.

        Parameters
        ----------
        mode : str
            The dust sublimation mode, which can be:
                * 'no'   - no sublimation
                * 'fast' - remove all dust in cells exceeding the
                           sublimation temperature
                * 'slow' - reduce the dust in cells exceeding the
                           sublimation temperature
                * 'cap'  - any temperature exceeding the sublimation
                           temperature is reset to the sublimation
                           temperature.

        temperature : float, optional
            The dust sublimation temperature, in K
        '''

        if mode not in ['no', 'fast', 'slow', 'cap']:
            raise Exception("mode should be one of no/fast/slow/cap")

        if mode != 'no' and temperature is None:
            raise Exception("Need to specify a sublimation temperature")

        self.sublimation_mode = mode
        self.sublimation_energy = self.temperature2specific_energy(temperature)

    def set_sublimation_specific_energy(self, mode, specific_energy=0.):
        '''
        Set the dust sublimation mode and specific energy.

        Parameters
        ----------
        mode : str
            The dust sublimation mode, which can be:
                * 'no'   - no sublimation
                * 'fast' - remove all dust in cells exceeding the
                           sublimation specific energy
                * 'slow' - reduce the dust in cells exceeding the
                           sublimation specific energy
                * 'cap'  - any specific energy exceeding the sublimation
                           specific energy is reset to the sublimation
                           specific energy.

        specific_energy : float, optional
            The dust sublimation specific energy, in cgs
        '''

        if mode not in ['no', 'fast', 'slow', 'cap']:
            raise Exception("mode should be one of no/fast/slow/cap")

        if mode != 'no' and specific_energy is None:
            raise Exception("Need to specify a sublimation specific_energy")

        self.sublimation_mode = mode
        self.sublimation_energy = specific_energy

    def _write_dust_sublimation(self, group):
        group.attrs['sublimation_mode'] = np.string_(self.sublimation_mode)
        if self.sublimation_mode in ['slow', 'fast', 'cap']:
            group.attrs['sublimation_specific_energy'] = self.sublimation_energy

    def _compute_mean_opacities(self):
        if not self.mean_opacities.all_set():
            self.mean_opacities.compute(self.optical_properties)

    def write(self, filename, compression=True):
        '''
        Write out to a standard dust file, including calculations of the mean
        opacities and optionally thermal emissivities.
        '''

        # Check that the optical properties have been set
        self.optical_properties.ensure_all_set()

        # Compute mean opacities if not already existent
        self._compute_mean_opacities()

        # Check that emissivities are set (before computing mean opacities)
        if not self.emissivities.all_set():
            logger.info("Computing emissivities assuming LTE")
            self.emissivities.set_lte(self.optical_properties, self.mean_opacities)

        # Create dust table set
        if isinstance(filename, six.string_types):
            dt = h5py.File(filename, 'w')
        else:
            dt = filename

        # Add standard keywords to header
        dt.attrs['version'] = 2
        dt.attrs['type'] = 1
        dt.attrs['python_version'] = np.string_(__version__)
        if self.md5:
            dt.attrs['asciimd5'] = np.string_(self.md5)

        # Add optical properties and scattering angle tables
        self.optical_properties.to_hdf5_group(dt)

        # Add mean opacities table
        self.mean_opacities.to_hdf5_group(dt)

        # Add emissivities and emissivity variable tables
        self.emissivities.to_hdf5_group(dt)

        # Dust sublimation parameters
        self._write_dust_sublimation(dt)

        # Close dust file
        if isinstance(dt, h5py.highlevel.File):
            dt.close()

        self._file = (filename, self.hash())

    def read(self, filename):
        '''
        Read in from a standard dust file
        '''

        from ..util.functions import asstr

        if isinstance(filename, six.string_types):

            # Check file exists
            if not os.path.exists(filename):
                raise Exception("File not found: %s" % filename)

            # Read in dust table set
            dt = h5py.File(filename, 'r')
            close = True

        else:

            # Read in dust table set
            dt = filename
            close = False

        # Check version and type
        if dt.attrs['version'] not in [1, 2]:
            raise Exception("Version should be 1 or 2")
        if dt.attrs['type'] != 1:
            raise Exception("Type should be 1")
        if 'asciimd5' in dt.attrs:
            self.md5 = asstr(dt.attrs['asciimd5'])
        else:
            self.md5 = None

        # Read in the optical properties
        self.optical_properties.from_hdf5_group(dt)

        # Read in the planck and rosseland mean opacities
        if dt.attrs['version'] == 1:
            logger.warn("Version 1 dust file detected - discarding mean opacities and recomputing them")
            self.mean_opacities.compute(self.optical_properties)
        else:
            self.mean_opacities.from_hdf5_group(dt)

        # Read in emissivities
        self.emissivities.from_hdf5_group(dt)

        # Close file object if needed
        if close:
            dt.close()
            self._file = (filename, self.hash())

    def chi_nu_temperature(self, temperature):
        """
        Compute the mean opacity to extinction for a blackbody at a given temperature.

        Parameters
        ----------
        temperature : float
            The temperature of the blackbody to use

        Returns
        -------
        chi_nu_mean : float
            The mean opacity to extinction
        """
        self._compute_mean_opacities()
        return interp1d_fast_loglog(self.mean_opacities.temperature,
                                    self.mean_opacities.chi_planck,
                                    temperature,
                                    bounds_error=True)

    def kappa_nu_temperature(self, temperature):
        """
        Compute the mean opacity to absorption for a blackbody at a given temperature.

        Parameters
        ----------
        temperature : float
            The temperature of the blackbody to use

        Returns
        -------
        kappa_nu_mean : float
            The mean opacity to absorption
        """
        self._compute_mean_opacities()
        return interp1d_fast_loglog(self.mean_opacities.temperature,
                                    self.mean_opacities.kappa_planck,
                                    temperature,
                                    bounds_error=True)

    def chi_nu_spectrum(self, nu, fnu):
        """
        Compute the mean opacity to extinction for a given spectrum.

        Parameters
        ----------
        nu : array_like
            The frequencies, in Hz
        fnu : array_like
            The monochromatic fluxes per unit frequency. Units are unimportant
            since proportionality constants are cancelled out in the
            computation.

        Returns
        -------
        chi_nu_mean : float
            The mean opacity to extinction
        """
        if nu.max() > self.optical_properties.nu.max() or nu.min() < self.optical_properties.nu.min():
            raise Exception("Opacity to extinction is not defined at all "
                            "spectrum frequencies")
        chi_nu = self.optical_properties.interp_chi_nu(nu)
        return (integrate_loglog(nu, fnu * chi_nu) /
                integrate_loglog(nu, fnu))

    def kappa_nu_spectrum(self, nu, fnu):
        """
        Compute the mean opacity to absorption for a given spectrum.

        Parameters
        ----------
        nu : array_like
            The frequencies, in Hz
        fnu : array_like
            The monochromatic fluxes per unit frequency. Units are unimportant
            since proportionality constants are cancelled out in the
            computation.

        Returns
        -------
        kappa_nu_mean : float
            The mean opacity to absorption
        """
        if nu.max() > self.optical_properties.nu.max() or nu.min() < self.optical_properties.nu.min():
            raise Exception("Opacity to absorption is not defined at all "
                            "spectrum frequencies")
        kappa_nu = self.optical_properties.interp_kappa_nu(nu)
        return (integrate_loglog(nu, fnu * kappa_nu) /
                integrate_loglog(nu, fnu))

    def temperature2specific_energy(self, temperature):
        """
        Convert a temperature to its corresponding specific energy value.

        Parameters
        ----------
        temperature : float or array_like
            The temperature to convert

        Returns
        -------
        specific_energy : float or array_like
            The specific energy corresponding to the input temperature
        """

        self._compute_mean_opacities()

        specific_energy = interp1d_fast_loglog(self.mean_opacities.temperature,
                                               self.mean_opacities.specific_energy,
                                               temperature,
                                               bounds_error=False,
                                               fill_value=np.nan)

        if np.isscalar(temperature):
            if temperature < self.mean_opacities.temperature[0]:
                specific_energy = self.mean_opacities.specific_energy[0]
            elif temperature > self.mean_opacities.temperature[-1]:
                specific_energy = self.mean_opacities.specific_energy[-1]
        else:
            specific_energy[temperature < self.mean_opacities.temperature[0]] = self.mean_opacities.specific_energy[0]
            specific_energy[temperature > self.mean_opacities.temperature[-1]] = self.mean_opacities.specific_energy[-1]

        return specific_energy

    def specific_energy2temperature(self, specific_energy):
        """
        Convert a specific energy value to its corresponding temperature.

        Parameters
        ----------
        specific_energy : float or array_like
            The specific energy to convert

        Returns
        -------
        temperature : float or array_like
            The temperature corresponding to the input specific energy
        """

        self._compute_mean_opacities()

        temperature = interp1d_fast_loglog(self.mean_opacities.specific_energy,
                                           self.mean_opacities.temperature,
                                           specific_energy,
                                           bounds_error=False,
                                           fill_value=np.nan)

        if np.isscalar(specific_energy):
            if specific_energy < self.mean_opacities.specific_energy[0]:
                temperature = self.mean_opacities.temperature[0]
            elif specific_energy > self.mean_opacities.specific_energy[-1]:
                temperature = self.mean_opacities.temperature[-1]
        else:
            temperature[specific_energy < self.mean_opacities.specific_energy[0]] = self.mean_opacities.temperature[0]
            temperature[specific_energy > self.mean_opacities.specific_energy[-1]] = self.mean_opacities.temperature[-1]

        return temperature


class IsotropicDust(SphericalDust):
    """
    This class should be used for dust properties that include isotropic
    scattering. The dust properties should be instatiated as::

        d = IsotropicDust(nu, albedo, chi)

    where ``nu``, ``albedo``, and ``chi`` are 1-D Numpy arrays containing the
    frequencies, albedo, and opacity to extinction respectively.
    """

    def __init__(self, nu, albedo, chi):

        SphericalDust.__init__(self)

        # Set cos(theta) grid for computing the scattering matrix elements
        self.optical_properties.mu = np.linspace(-1., 1., 2)

        # Set optical properties
        self.optical_properties.nu = nu
        self.optical_properties.albedo = albedo
        self.optical_properties.chi = chi

        # Compute scattering matrix elements
        self.optical_properties.initialize_scattering_matrix()

        # Set scattering matrix to isotropic values
        self.optical_properties.P1[:, :] = 1.
        self.optical_properties.P2[:, :] = 0.
        self.optical_properties.P3[:, :] = 1.
        self.optical_properties.P4[:, :] = 0.

        # Sort optical properties
        self.optical_properties._sort()


class HenyeyGreensteinDust(SphericalDust):
    """
    This class should be used for dust properties that include
    scattering parameterized by the `Henyey-Greenstein, 1941
    <http://dx.doi.org/10.1086/144246>`_ function. The dust properties should
    be instatiated as::

        d = HenyeyGreensteinDust(nu, albedo, chi, g, p_lin_max)

    where ``nu``, ``albedo``, and ``chi`` are 1-D Numpy arrays containing the
    frequencies, albedo, and opacity to extinction respectively, and ``g`` and
    ``p_lin_max`` are also 1-D Numpy arrays containing the asymmetry parameter
    and the maximum linear polarization.
    """


    def __init__(self, nu, albedo, chi, g, p_lin_max):

        SphericalDust.__init__(self)

        # Set cos(theta) grid for computing the scattering matrix elements
        n_mu = 100
        self.optical_properties.mu = np.linspace(-1., 1., n_mu)

        # Set optical properties
        self.optical_properties.nu = nu
        self.optical_properties.albedo = albedo
        self.optical_properties.chi = chi

        # Compute scattering matrix elements
        self.optical_properties.initialize_scattering_matrix()

        for i in range(n_mu):
            self.optical_properties.P1[:, i], \
            self.optical_properties.P2[:, i], \
            self.optical_properties.P3[:, i], \
            self.optical_properties.P4[:, i] = henyey_greenstein(self.optical_properties.mu[i], g, p_lin_max)


class HOCHUNKDust(HenyeyGreensteinDust):
    """
    This class should be used for dust properties that include
    scattering parameterized by the `Henyey-Greenstein, 1941
    <http://dx.doi.org/10.1086/144246>`_ function, which are formatted for the
    `HOCHUNK code <http://gemelli.colorado.edu/~bwhitney/codes/>`_. The dust
    properties should be instatiated as::

        d = HOCHUNKDust(filename)

    where ``filename`` is the name of the file containing the dust properties
    in the HOCHUNK format.
    """
    def __init__(self, filename):

        # Read in dust file
        dustfile = np.loadtxt(filename, dtype=[('wav', float), ('c_ext', float), \
                              ('c_sca', float), ('chi', float), ('g', float), \
                              ('p_lin_max', float)], usecols=[0, 1, 2, 3, 4, 5])

        # Ensure file is ordered in increasing frequency
        if dustfile['wav'][-1] > dustfile['wav'][0]:
            dustfile = dustfile[::-1]

        # Compute frequency and albedo
        nu = c / dustfile['wav'] * 1.e4
        albedo = dustfile['c_sca'] / dustfile['c_ext']

        self.md5 = hashlib.md5(open(filename, 'rb').read()).hexdigest()

        HenyeyGreensteinDust.__init__(self, nu, albedo, dustfile['chi'], dustfile['g'], dustfile['p_lin_max'])

TTsreDust = HOCHUNKDust

class CoatsphSingle(SphericalDust):

    def __init__(self, directory, size, density):
        '''
        Initialize single-component dust.

        Parameters
        ----------
        directory : str
            Directory containing all the files describing the dust
        size : float
            Grain size, in cm
        density : float
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

        n_wav = len(dustfile)

        self.optical_properties.nu = c / dustfile['wav'] * 1.e4
        self.optical_properties.albedo = dustfile['q_sca'] / dustfile['q_ext']
        self.optical_properties.chi = 0.75 * dustfile['q_ext'] / size / density

        # Read in scattering matrix elements

        for i in range(n_wav):

            filename = '%s/coatsph_scat_%04i_0001.dat' % (directory, i + 1)

            phasefile = np.loadtxt(filename, skiprows=9,
                        dtype=[('theta', float), ('s11', float), ('polariz',
                        float), ('s12', float), ('s33', float), ('s34',
                        float)])

            if i == 0:
                self.optical_properties.mu = np.cos(np.radians(phasefile['theta']))
                self.optical_properties.initialize_scattering_matrix()

            self.optical_properties.P1[i, :] = phasefile['s11']
            self.optical_properties.P2[i, :] = phasefile['s12']
            self.optical_properties.P3[i, :] = phasefile['s33']
            self.optical_properties.P4[i, :] = phasefile['s34']


class CoatsphMultiple(SphericalDust):

    def __init__(self, directory):
        '''
        Initialize multi-component dust.

        Parameters
        ----------
        directory : str
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

        n_wav = len(dustfile)
        self.optical_properties.nu = c / dustfile['wav'] * 1.e4
        self.optical_properties.albedo = dustfile['c_sca'] / dustfile['c_ext']
        self.optical_properties.chi = dustfile['chi']

        # Read in scattering matrix elements

        for i in range(n_wav):

            filename = '%s/coatsph_scat.%04i.dat' % (directory, i + 1)

            phasefile = np.loadtxt(filename, skiprows=7,
                        dtype=[('theta', float), ('s11', float), ('polariz',
                        float), ('s12', float), ('s33', float), ('s34',
                        float)])

            if i == 0:
                self.optical_properties.mu = np.cos(np.radians(phasefile['theta']))
                self.optical_properties.initialize_scattering_matrix()

            self.optical_properties.P1[i, :] = phasefile['s11']
            self.optical_properties.P2[i, :] = phasefile['s12']
            self.optical_properties.P3[i, :] = phasefile['s33']
            self.optical_properties.P4[i, :] = phasefile['s34']


class MieXDust(SphericalDust):

    def __init__(self, model):

        SphericalDust.__init__(self)

        wav = np.loadtxt('%s.alb' % model, usecols=[0])
        self.optical_properties.albedo = np.loadtxt('%s.alb' % model, usecols=[1])
        kappa = np.loadtxt('%s.k_abs' % model, usecols=[1])
        self.optical_properties.chi = kappa / (1 - self.optical_properties.albedo)

        # Check for NaN values
        for quantity in ['chi', 'albedo']:

            values = self.optical_properties.__dict__[quantity]

            if np.any(np.isnan(values)):
                logger.warn("NaN values found inside MieX %s file - interpolating" % quantity)
                invalid = np.isnan(values)
                values[invalid] = interp1d_fast_loglog(wav[~invalid], values[~invalid], wav[invalid])
                if np.any(np.isnan(values)):
                    raise Exception("Did not manage to fix NaN values in MieX %s" % quantity)

        self.optical_properties.nu = c / wav * 1.e4

        n_wav = len(wav)
        n_mu = (len(open('%s.f11' % model).readlines()) // n_wav) - 1

        mu = np.zeros(n_mu)

        # Read mu
        f11 = open('%s.f11' % model)
        f11.readline()
        f11.readline()
        for i in range(n_mu):
            mu[i] = np.cos(np.radians(float(f11.readline().split()[0])))
        f11.close()
        self.optical_properties.mu = mu[::-1]

        # Read in matrix elements

        self.optical_properties.initialize_scattering_matrix()

        f11 = open('%s.f11' % model)
        f12 = open('%s.f12' % model)
        f33 = open('%s.f33' % model)
        f34 = open('%s.f34' % model)

        f11.readline()
        f12.readline()
        f33.readline()
        f34.readline()

        for j in range(n_wav):

            if float(f11.readline()) != wav[j]:
                raise Exception("Incorrect wavelength in f11")
            if float(f12.readline()) != wav[j]:
                raise Exception("Incorrect wavelength in f12")
            if float(f33.readline()) != wav[j]:
                raise Exception("Incorrect wavelength in f33")
            if float(f34.readline()) != wav[j]:
                raise Exception("Incorrect wavelength in f34")

            for i in range(n_mu):

                self.optical_properties.P1[j, n_mu - i - 1] = float(f11.readline().split()[1])
                self.optical_properties.P2[j, n_mu - i - 1] = float(f12.readline().split()[1])
                self.optical_properties.P3[j, n_mu - i - 1] = float(f33.readline().split()[1])
                self.optical_properties.P4[j, n_mu - i - 1] = float(f34.readline().split()[1])

        for i in range(n_mu):

            for quantity in ['P1', 'P2', 'P3', 'P4']:

                values = self.optical_properties.__dict__[quantity]

                if np.any(np.isnan(values[:, i])):
                    logger.warn("NaN values found inside MieX %s file - interpolating" % quantity)
                    invalid = np.isnan(values[:, i])
                    values[:, i][invalid] = interp1d_fast_loglog(wav[~invalid], values[:, i][~invalid], wav[invalid])
                    if np.any(np.isnan(values[:, i])):
                        raise Exception("Did not manage to fix NaN values in MieX %s" % quantity)


class BHDust(SphericalDust):
    """
    This class should be used for dust properties that were computed using
    `this dust calculation code <https://github.com/hyperion-rt/bhmie>`_ which
    is a wrapper to the ``bhmie`` routine originally written by C.F. Bohren and
    D. Huffman and improved by B. Draine.

    When using the ``bhmie`` code, you should set the output format to ``2``,
    which will create a number of files ending in ``.wav``, ``.mu``, ``.alb``,
    etc. Then, instantiate this class with the name of the directory containing
    these output files along with the prefix used. For example, if you use
    ``directory/mydust`` as a prefix in ``bhmie``, you can import this dust
    with::

        >>> from hyperion.dust import BHDust
        >>> d = BHDust('directory/mydust')

    """
    def __init__(self, model):

        SphericalDust.__init__(self)

        mu = np.loadtxt('%s.mu' % model)

        nu = c / np.loadtxt('%s.wav' % model) * 1.e4
        albedo = np.loadtxt('%s.alb' % model)
        chi = np.loadtxt('%s.chi' % model)

        P1 = np.loadtxt('%s.f11' % model)
        P2 = np.loadtxt('%s.f12' % model)
        P3 = np.loadtxt('%s.f33' % model)
        P4 = np.loadtxt('%s.f34' % model)

        if nu[-1] < nu[0]:
            nu = nu[::-1]
            albedo = albedo[::-1]
            chi = chi[::-1]
            P1 = P1[::-1, :]
            P2 = P2[::-1, :]
            P3 = P3[::-1, :]
            P4 = P4[::-1, :]

        if mu[-1] < mu[0]:
            mu = mu[::-1]
            P1 = P1[:, ::-1]
            P2 = P2[:, ::-1]
            P3 = P3[:, ::-1]
            P4 = P4[:, ::-1]

        self.optical_properties.mu = mu

        self.optical_properties.nu = nu
        self.optical_properties.albedo = albedo
        self.optical_properties.chi = chi

        self.optical_properties.P1 = P1
        self.optical_properties.P2 = P2
        self.optical_properties.P3 = P3
        self.optical_properties.P4 = P4
