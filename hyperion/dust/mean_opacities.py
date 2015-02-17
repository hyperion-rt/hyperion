from __future__ import print_function, division

import hashlib

import numpy as np
from astropy.table import Table, Column
from astropy import log as logger

from ..util.integrate import integrate_loglog
from ..util.interpolate import interp1d_fast_loglog
from ..util.functions import FreezableClass, nu_common, B_nu, dB_nu_dT, planck_nu_range
from ..util.constants import sigma


class MeanOpacities(FreezableClass):

    def __init__(self):

        self.specific_energy = None
        self.temperature = None
        self.chi_planck = None
        self.kappa_planck = None
        self.chi_inv_planck = None
        self.kappa_inv_planck = None
        self.chi_rosseland = None
        self.kappa_rosseland = None
        self._freeze()

    def compute(self, optical_properties, n_temp=1200, temp_min=0.1, temp_max=100000.):
        """
        Compute various mean opacities:

            * Planck mean opacity
            * Reciprocal Planck mean opacity
            * Rosseland mean opacity
        """

        # Set temperatures to compute the mean opacities for
        temperatures = np.logspace(np.log10(temp_min),
                                   np.log10(temp_max), n_temp)

        # To avoid issues that may be confusing to users if they ask for
        # temperatures at exactly the min max, we reset the temperature min/max
        # manually (otherwise the np.log10 and subsequent 10** cause a loss in
        # precision)
        temperatures[0] = temp_min
        temperatures[-1] = temp_max

        # Find common frequency scale
        planck_nu = planck_nu_range(temp_min, temp_max)
        nu = nu_common(planck_nu, optical_properties.nu)


        if planck_nu.min() < optical_properties.nu.min():
            logger.warn("Planck function for lowest temperature not completely covered by opacity function")
            nu = nu[nu >= optical_properties.nu.min()]

        if planck_nu.max() > optical_properties.nu.max():
            logger.warn("Planck function for highest temperature not completely covered by opacity function")
            nu = nu[nu <= optical_properties.nu.max()]

        # Interpolate opacity to new frequency grid
        chi_nu = interp1d_fast_loglog(optical_properties.nu,
                                      optical_properties.chi, nu)
        kappa_nu = interp1d_fast_loglog(optical_properties.nu,
                                        optical_properties.kappa, nu)

        # Set mean opacity variable
        self.var_name = 'specific_energy'

        # Initialize mean opacity arrays
        self.chi_planck = np.zeros(n_temp)
        self.kappa_planck = np.zeros(n_temp)
        self.chi_inv_planck = np.zeros(n_temp)
        self.kappa_inv_planck = np.zeros(n_temp)
        self.chi_rosseland = np.zeros(n_temp)
        self.kappa_rosseland = np.zeros(n_temp)

        # Loop through the emissivities and compute mean opacities
        for it, T in enumerate(temperatures):

            # Compute Planck function and derivative with respect to temperature
            b_nu = B_nu(nu, T)
            db_nu_dt = dB_nu_dT(nu, T)

            # Compute planck mean opacity
            self.chi_planck[it] = (integrate_loglog(nu, b_nu * chi_nu) /
                                   integrate_loglog(nu, b_nu))

            # Compute planck mean absoptive opacity
            self.kappa_planck[it] = (integrate_loglog(nu, b_nu * kappa_nu) /
                                     integrate_loglog(nu, b_nu))

            # Compute reciprocal planck mean opacity
            self.chi_inv_planck[it] = (integrate_loglog(nu, b_nu) /
                                       integrate_loglog(nu, b_nu / chi_nu))

            # Compute reciprocal planck mean aborptive opacity
            self.kappa_inv_planck[it] = (integrate_loglog(nu, b_nu) /
                                         integrate_loglog(nu, b_nu / kappa_nu))

            # Compute rosseland mean opacity
            self.chi_rosseland[it] = (integrate_loglog(nu, db_nu_dt) /
                                      integrate_loglog(nu, db_nu_dt / chi_nu))

            # Compute rosseland mean aborptive opacity
            self.kappa_rosseland[it] = (integrate_loglog(nu, db_nu_dt) /
                                        integrate_loglog(nu, db_nu_dt / kappa_nu))

        self.temperature = temperatures
        self.specific_energy = 4. * sigma * temperatures ** 4. * self.kappa_planck

    def to_hdf5_group(self, group):

        if not self.all_set():
            raise Exception("Not all attributes of the mean opacities are set")

        # Create mean opacities table
        tmean = Table()
        tmean.add_column(Column(data=self.temperature, name='temperature'))
        tmean.add_column(Column(data=self.specific_energy, name='specific_energy'))
        tmean.add_column(Column(data=self.chi_planck, name='chi_planck'))
        tmean.add_column(Column(data=self.kappa_planck, name='kappa_planck'))
        tmean.add_column(Column(data=self.chi_inv_planck, name='chi_inv_planck'))
        tmean.add_column(Column(data=self.kappa_inv_planck, name='kappa_inv_planck'))
        tmean.add_column(Column(data=self.chi_rosseland, name='chi_rosseland'))
        tmean.add_column(Column(data=self.kappa_rosseland, name='kappa_rosseland'))

        # Add to group
        tmean.write(group, path='mean_opacities')

    def from_hdf5_group(self, group):

        tmean = Table.read(group, path='mean_opacities')
        self.temperature = tmean['temperature']
        self.specific_energy = tmean['specific_energy']
        self.chi_planck = tmean['chi_planck']
        self.kappa_planck = tmean['kappa_planck']
        self.chi_inv_planck = tmean['chi_inv_planck']
        self.kappa_inv_planck = tmean['kappa_inv_planck']
        self.chi_rosseland = tmean['chi_rosseland']
        self.kappa_rosseland = tmean['kappa_rosseland']

    def all_set(self):
        return (self.temperature is not None and
                self.specific_energy is not None and
                self.chi_planck is not None and
                self.kappa_planck is not None and
                self.chi_inv_planck is not None and
                self.kappa_inv_planck is not None and
                self.chi_rosseland is not None and
                self.kappa_rosseland is not None)

    def plot(self, figure, subplot):

        if not self.all_set():
            raise Exception("Not all attributes of the mean opacities are set")

        ax = figure.add_subplot(subplot)

        ax.loglog(self.specific_energy, self.chi_planck, color='red', label='Planck Extinction')
        ax.loglog(self.specific_energy, self.kappa_planck, color='orange', label='Planck Absorption')
        ax.loglog(self.specific_energy, self.chi_inv_planck, color='blue', label='Reciprocal Planck Extinction')
        ax.loglog(self.specific_energy, self.kappa_inv_planck, color='lightblue', label='Reciprocal Planck Absorption')
        ax.loglog(self.specific_energy, self.chi_rosseland, color='green', label='Rosseland Extinction')
        ax.loglog(self.specific_energy, self.kappa_rosseland, color='lightgreen', label='Rosseland Absorption')
        ax.legend(loc=2)
        ax.set_xlabel("Specific energy (ergs/s/g)")
        ax.set_ylabel("Mean opacity (cm^2/g)")
        ax.set_xlim(self.specific_energy.min(), self.specific_energy.max())

        return figure

    def hash(self):
        h = hashlib.md5()
        h.update(self.temperature.tostring())
        h.update(self.specific_energy.tostring())
        h.update(self.chi_planck.tostring())
        h.update(self.kappa_planck.tostring())
        h.update(self.chi_inv_planck.tostring())
        h.update(self.kappa_inv_planck.tostring())
        h.update(self.chi_rosseland.tostring())
        h.update(self.kappa_rosseland.tostring())
        return h.hexdigest()
