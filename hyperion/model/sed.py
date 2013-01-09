import numpy as np

from ..util.functions import FreezableClass
from ..util.constants import c


class SED(FreezableClass):
    """
    Class to represent an SED or set of SEDs

    Parameters
    ----------
    nu : ndarray
        The frequencies at which the SED is defined, in Hz
    flux : ndarray, optional
        The fluxes for the SED. The last dimensions should match the number of
        frequencies. This flux can be f_nu or nu * f_nu.
    unc : ndarray, optional
        The flux uncertainties for the SED. The last dimensions should match
        the number of frequencies.
    units : str
        The units of the flux
    """

    def __init__(self, nu, flux=None, unc=None, units=None):
        self.nu = nu
        self.flux = flux
        self.unc = unc
        self.units = units
        self._freeze()

    @property
    def nu(self):
        return self._nu

    @nu.setter
    def nu(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._nu = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            self._nu = value
        else:
            raise TypeError("nu should be a 1-d sequence")

    @property
    def flux(self):
        return self._flux

    @flux.setter
    def flux(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._flux = value
        elif isinstance(value, np.ndarray) and value.ndim >= 1:
            if self.nu is not None and len(self.nu) != value.shape[-1]:
                raise ValueError("the last dimension of the flux array should match the length of the nu array (expected {0} but found {1})".format(len(self.nu), value.shape[-1]))
            else:
                if hasattr(self, 'unc') and self.unc is not None:
                    if value.shape != self.unc.shape:
                        raise ValueError("dimensions should match that of unc")
                self._flux = value
        else:
            raise TypeError("flux should be a multi-dimensional array")

    @property
    def unc(self):
        return self._unc

    @unc.setter
    def unc(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._unc = value
        elif isinstance(value, np.ndarray) and value.ndim >= 1:
            if self.nu is not None and len(self.nu) != value.shape[-1]:
                raise ValueError("the last dimension of the unc array should match the length of the nu array (expected {0} but found {1})".format(len(self.nu), value.shape[-1]))
            else:
                if hasattr(self, 'flux') and  self.flux is not None:
                    if value.shape != self.flux.shape:
                        raise ValueError("dimensions should match that of flux")
                self._unc = value
        else:
            raise TypeError("unc should be a multi-dimensional array")
    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if value is None or isinstance(value, basestring):
            self._unit = value
        else:
            raise ValueError("unit should be a string")

    @property
    def wav(self):
        return c / self.nu * 1e4

    def __iter__(self):
        if self.unc is None:
            return (x for x in [self.wav, self.flux])
        else:
            return (x for x in [self.wav, self.flux, self.unc])
