from __future__ import print_function, division

import hashlib

import numpy as np
from astropy.table import Table, Column

from ..util.integrate import integrate_loglog
from ..util.interpolate import interp1d_fast_loglog
from ..util.functions import FreezableClass, nu_common
from ..util.constants import sigma


class MeanOpacities(FreezableClass):

    def __init__(self):

        self.var_name = None
        self.var = None
        self.chi_planck = None
        self.kappa_planck = None
        self.chi_rosseland = None
        self.kappa_rosseland = None
        self._freeze()

    def compute(self, emissivities, optical_properties):

        # Find common frequency scale
        nu = nu_common(emissivities.nu, optical_properties.nu)

        # Interpolate opacity to new frequency grid
        chi_nu = interp1d_fast_loglog(optical_properties.nu,
                                      optical_properties.chi, nu)
        kappa_nu = interp1d_fast_loglog(optical_properties.nu,
                                        optical_properties.kappa, nu)

        # Find number of emissivities to compute mean opacities for
        n_emiss = len(emissivities.var)

        # Set mean opacity variable
        self.var_name = 'specific_energy'
        self.var = emissivities.var

        # Initialize mean opacity arrays
        self.chi_planck = np.zeros(n_emiss)
        self.kappa_planck = np.zeros(n_emiss)
        self.chi_rosseland = np.zeros(n_emiss)
        self.kappa_rosseland = np.zeros(n_emiss)

        # Loop through the emissivities and compute mean opacities
        for ivar in range(n_emiss):

            # Extract emissivity and interpolate to common frequency array
            jnu = interp1d_fast_loglog(emissivities.nu,
                                       emissivities.jnu[:, ivar], nu,
                                       bounds_error=False, fill_value=0.)

            # Define I_nu = J_nu / kappa_nu
            inu = jnu / kappa_nu

            # Compute planck mean opacity
            self.chi_planck[ivar] = integrate_loglog(nu, inu * chi_nu) \
                                  / integrate_loglog(nu, inu)

            # Compute planck mean absoptive opacity
            self.kappa_planck[ivar] = integrate_loglog(nu, inu * kappa_nu) \
                                    / integrate_loglog(nu, inu)

            # Compute Rosseland mean opacity
            self.chi_rosseland[ivar] = integrate_loglog(nu, inu) \
                                     / integrate_loglog(nu, inu / chi_nu)

            # Compute Rosseland mean opacity
            self.kappa_rosseland[ivar] = integrate_loglog(nu, inu) \
                                       / integrate_loglog(nu, inu / kappa_nu)

    def _temperature2specific_energy(self, temperature):
        temperatures = np.sqrt(np.sqrt((self.var / (4. * sigma * self.kappa_planck))))
        specific_energy = interp1d_fast_loglog(temperatures, self.var, temperature,
                                               bounds_error=False, fill_value=np.nan)
        if np.isscalar(temperature):
            if temperature < temperatures[0]:
                specific_energy = self.var[0]
            elif temperature > temperatures[-1]:
                specific_energy = self.var[-1]
        else:
            specific_energy[temperature < temperatures[0]] = self.var[0]
            specific_energy[temperature > temperatures[-1]] = self.var[-1]
        return specific_energy

    def _specific_energy2temperature(self, specific_energy):
        temperatures = np.sqrt(np.sqrt((self.var / (4. * sigma * self.kappa_planck))))
        temperature = interp1d_fast_loglog(self.var, temperatures, specific_energy,
                                           bounds_error=False, fill_value=np.nan)
        if np.isscalar(specific_energy):
            if specific_energy < self.var[0]:
                temperature = temperatures[0]
            elif specific_energy > self.var[-1]:
                temperature = temperatures[-1]
        else:
            temperature[specific_energy < self.var[0]] = temperatures[0]
            temperature[specific_energy > self.var[-1]] = temperatures[-1]
        return temperature

    def to_hdf5_group(self, group):

        if not self.all_set():
            raise Exception("Not all attributes of the mean opacities are set")

        # Create mean opacities table
        tmean = Table()
        tmean.meta['var_name'] = np.string_(self.var_name)
        tmean.add_column(Column(data=self.var, name=self.var_name,))
        tmean.add_column(Column(data=self.chi_planck, name='chi_planck'))
        tmean.add_column(Column(data=self.kappa_planck, name='kappa_planck'))
        tmean.add_column(Column(data=self.chi_rosseland, name='chi_rosseland'))
        tmean.add_column(Column(data=self.kappa_rosseland, name='kappa_rosseland'))

        # Add to group
        tmean.write(group, path='mean_opacities')

    def from_hdf5_group(self, group):

        from ..util.functions import asstr

        tmean = Table.read(group, path='mean_opacities')
        self.var_name = asstr(tmean.meta['var_name'])
        self.var = tmean[self.var_name]
        self.chi_planck = tmean['chi_planck']
        self.kappa_planck = tmean['kappa_planck']
        self.chi_rosseland = tmean['chi_rosseland']
        self.kappa_rosseland = tmean['kappa_rosseland']

    def all_set(self):
        return self.var_name is not None and \
               self.var is not None and \
               self.chi_planck is not None and \
               self.kappa_planck is not None and \
               self.chi_rosseland is not None and \
               self.kappa_rosseland is not None

    def plot(self, figure, subplot):

        if not self.all_set():
            raise Exception("Not all attributes of the mean opacities are set")

        ax = figure.add_subplot(subplot)

        ax.loglog(self.var, self.chi_planck, color='red', label='Planck Extinction')
        ax.loglog(self.var, self.kappa_planck, color='orange', label='Planck Absorption')
        ax.loglog(self.var, self.chi_rosseland, color='blue', label='Rosseland Extinction')
        ax.loglog(self.var, self.kappa_rosseland, color='lightblue', label='Rosseland Absorption')
        ax.legend(loc=2)
        ax.set_xlabel("Specific energy (ergs/s/g)")
        ax.set_ylabel("Mean opacity (cm^2/g)")
        ax.set_xlim(self.var.min(), self.var.max())

        return figure

    def hash(self):
        h = hashlib.md5()
        h.update(self.var_name.encode('utf-8'))
        h.update(self.var.tostring())
        h.update(self.chi_planck.tostring())
        h.update(self.kappa_planck.tostring())
        h.update(self.chi_rosseland.tostring())
        h.update(self.kappa_rosseland.tostring())
        return h.hexdigest()
