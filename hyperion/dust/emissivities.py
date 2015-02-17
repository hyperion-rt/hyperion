from __future__ import print_function, division

import hashlib

import numpy as np
from astropy.table import Table, Column

from ..util.integrate import integrate_loglog
from ..util.interpolate import interp1d_fast_loglog
from ..util.functions import B_nu, FreezableClass, nu_common, \
                                    planck_nu_range, bool2str, is_numpy_array, monotonically_increasing
from astropy import log as logger


class Emissivities(FreezableClass):

    def __init__(self):

        self.is_lte = False
        self.var_name = None
        self.var = None
        self.nu = None
        self.jnu = None

        self._freeze()

    def normalize(self):
        for ivar in range(len(self.var)):
            norm = integrate_loglog(self.nu, self.jnu[:, ivar] / self.nu)
            self.jnu[:, ivar] /= norm

    def set_lte(self, optical_properties, mean_opacities):

        # Specify that emissivities are LTE
        self.is_lte = True

        # Get temperatures from mean opacities
        temperature = mean_opacities.temperature
        specific_energy = mean_opacities.specific_energy

        # Set frequency scale
        planck_nu = planck_nu_range(temperature[0], temperature[-1])
        self.nu = nu_common(planck_nu, optical_properties.nu)

        if planck_nu.min() < optical_properties.nu.min():
            logger.warn("Planck function for lowest temperature not completely covered by opacity function")
            self.nu = self.nu[self.nu >= optical_properties.nu.min()]

        if planck_nu.max() > optical_properties.nu.max():
            logger.warn("Planck function for highest temperature not completely covered by opacity function")
            self.nu = self.nu[self.nu <= optical_properties.nu.max()]

        # Compute opacity to absorption
        kappa_nu = interp1d_fast_loglog(optical_properties.nu,
                                        optical_properties.kappa, self.nu)

        # Compute LTE emissivities
        self.var_name = 'specific_energy'
        self.var = specific_energy
        self.jnu = np.zeros((len(self.nu), len(temperature)))

        # Find LTE emissivities
        for it, T in enumerate(temperature):
            self.jnu[:, it] = kappa_nu * B_nu(self.nu, T)

    def to_hdf5_group(self, group):

        if not self.all_set():
            raise Exception("Not all attributes of the emissivities are set")

        # Write out the emissivity variable type
        if self.var_name == 'specific_energy':
            group.attrs['emissvar'] = np.string_('E')
        else:
            raise Exception("Unknown emissivity variable: %s" % self.var_name)

        # Create emissivity variable table
        temissvar = Table()
        temissvar.add_column(Column(data=self.var, name=self.var_name))

        # Create emissivities table
        temiss = Table()
        temiss.add_column(Column(data=self.nu, name='nu'))
        temiss.add_column(Column(data=self.jnu, name='jnu'))

        group.attrs['lte'] = bool2str(self.is_lte)

        # Add to group
        temissvar.write(group, path='emissivity_variable')
        temiss.write(group, path='emissivities')

    def from_hdf5_group(self, group):

        from ..util.functions import asstr

        # Find the emissivity variable type
        if asstr(group.attrs['emissvar']) == 'E':
            self.var_name = 'specific_energy'
        else:
            raise Exception("Unknown emissivity variable: %s" %
                            group.attrs['emissvar'])

        # Read in emissivity variable
        temissvar = Table.read(group, path='emissivity_variable')
        self.var = temissvar[self.var_name]

        # Read emissivities
        temiss = Table.read(group, path='emissivities')
        self.nu = temiss['nu']
        self.jnu = temiss['jnu']
        self.is_lte = group.attrs['lte'].decode('utf-8').lower() == 'yes'

    def all_set(self):
        return self.var_name is not None and \
               self.var is not None and \
               self.nu is not None and \
               self.jnu is not None

    def plot(self, figure, subplot):

        if not self.all_set():
            raise Exception("Not all attributes of the emissivities are set")

        import matplotlib.pyplot as plt

        self.normalize()
        peak = self.jnu.max()

        m = plt.cm.gist_heat
        vmin, vmax = np.log10(peak) - 6., np.log10(peak)

        ax = figure.add_subplot(subplot)
        ax.patch.set_facecolor('black')
        ax.contourf(self.nu, self.var,
                     np.log10(np.clip(np.abs(self.jnu.swapaxes(0, 1)), 10. ** vmin, 10. ** vmax)),
                     np.linspace(vmin, vmax, 30),
                     cmap=m)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(self.nu.min(), self.nu.max())
        ax.set_ylim(self.var.min(), self.var.max())
        ax.set_title('Emissivities', y=0.9, verticalalignment='top',
                     color='white')

        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Specific energy (ergs/s/g)")

        return figure

    def hash(self):
        h = hashlib.md5()
        h.update(str(self.is_lte).encode('utf-8'))
        h.update(self.var_name.encode('utf-8'))
        h.update(self.var.tostring())
        h.update(self.nu.tostring())
        h.update(self.jnu.tostring())
        return h.hexdigest()

    def __setattr__(self, attribute, value):
        if attribute in ['nu', 'var'] and value is not None:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError(attribute + " should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError(attribute + " should be monotonically increasing")
            if value[0] <= 0.:
                raise ValueError(attribute + ' should be strictly positive')
        if attribute in ['jnu'] and value is not None:
            if self.nu is None:
                raise ValueError("nu needs to be set before " + attribute)
            if self.var is None:
                raise ValueError("var needs to be set before " + attribute)
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 2:
                raise ValueError(attribute + " should be a 2-D array")
            if value.shape[0] != len(self.nu) or value.shape[1] != len(self.var):
                raise ValueError(attribute + " has an incorrect shape: %s but expected (%i, %i)" % (value.shape, len(self.nu), len(self.var)))
            if np.any(value < 0.):
                raise ValueError("jnu should be positive")
        FreezableClass.__setattr__(self, attribute, value)
