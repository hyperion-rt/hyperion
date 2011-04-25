import atpy
import numpy as np

from hyperion.util.integrate import integrate_loglog
from hyperion.util.interpolate import interp1d_fast_loglog
from hyperion.util.functions import FreezableClass, nu_common


class MeanOpacities(FreezableClass):

    def __init__(self):

        self.set = False
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
        self.var_name = 'specific_energy_abs'
        self.var = emissivities.var

        # Initialize mean opacity arrays
        self.chi_planck = np.zeros(n_emiss)
        self.kappa_planck = np.zeros(n_emiss)
        self.chi_rosseland = np.zeros(n_emiss)
        self.kappa_rosseland = np.zeros(n_emiss)

        # Loop through the temperatures and compute mean opacities
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

        # Indicate that mean opacities have been set
        self.set = True

    def to_table_set(self, table_set):

        # Create mean opacities table
        tmean = atpy.Table(name='Mean opacities')
        tmean.add_keyword('var_name', self.var_name)
        tmean.add_column(self.var_name, self.var)
        tmean.add_column('chi_planck', self.chi_planck)
        tmean.add_column('kappa_planck', self.kappa_planck)
        tmean.add_column('chi_rosseland', self.chi_rosseland)
        tmean.add_column('kappa_rosseland', self.kappa_rosseland)

        # Add to table set
        table_set.append(tmean)

    def from_table_set(self, table_set):

        tmean = table_set['Mean opacities']
        self.var_name = tmean.keywords['var_name']
        self.specific_energy_abs = tmean[self.var_name]
        self.chi_planck = tmean['chi_planck']
        self.kappa_planck = tmean['kappa_planck']
        self.chi_rosseland = tmean['chi_rosseland']
        self.kappa_rosseland = tmean['kappa_rosseland']

    def plot(self, figure, subplot):

        ax = figure.add_subplot(subplot)

        ax.loglog(self.var, self.chi_planck, color='red')
        ax.loglog(self.var, self.kappa_planck, color='orange')
        ax.loglog(self.var, self.chi_rosseland, color='blue')
        ax.loglog(self.var, self.kappa_rosseland, color='lightblue')
        ax.set_xlim(self.var.min(), self.var.max())

        return figure
