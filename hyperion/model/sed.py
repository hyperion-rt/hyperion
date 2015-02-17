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
    val : ndarray, optional
        The values for the SED. The last dimensions should match the number of
        frequencies.
    unc : ndarray, optional
        The uncertainties for the SED values. The last dimensions should match
        the number of frequencies.
    units : str
        The units of the values
    """

    def __init__(self, nu, val=None, unc=None, units=None):

        self.nu = nu
        self.val = val
        self.unc = unc
        self.units = units

        self.ap_min = None
        self.ap_max = None

        self.d_min = None
        self.d_max = None

        self.distance = None

        self.inside_observer = False

        self._freeze()

    @property
    def nu(self):
        """
        The frequencies for which the SED is defined (in Hz)
        """
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
    def val(self):
        """
        The SED values (fluxes, flux densities, surface brightness, or polarization) in the units given by the ``.unit`` property.
        """
        return self._val

    @val.setter
    def val(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._val = value
        elif isinstance(value, np.ndarray) and value.ndim >= 1:
            if self.nu is not None and len(self.nu) != value.shape[-1]:
                raise ValueError("the last dimension of the value array should match the length of the nu array (expected {0} but found {1})".format(len(self.nu), value.shape[-1]))
            else:
                if hasattr(self, 'unc') and self.unc is not None:
                    if value.shape != self.unc.shape:
                        raise ValueError("dimensions should match that of unc")
                self._val = value
        else:
            raise TypeError("val should be a multi-dimensional array")

    @property
    def unc(self):
        """
        The uncertainties on the SED values in the units given by the ``.unit`` property.
        """
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
                if hasattr(self, 'val') and  self.val is not None:
                    if value.shape != self.val.shape:
                        raise ValueError("dimensions should match that of val")
                self._unc = value
        else:
            raise TypeError("unc should be a multi-dimensional array")

    @property
    def unit(self):
        """
        The units of the SED values.
        """
        return self._unit

    @unit.setter
    def unit(self, value):
        if value is None or isinstance(value, basestring):
            self._unit = value
        else:
            raise ValueError("unit should be a string")

    @property
    def wav(self):
        """
        The wavelengths for which the SED is defined (in microns).
        """
        return c / self.nu * 1e4

    def __iter__(self):
        if self.unc is None:
            return (x for x in [self.wav, self.val])
        else:
            return (x for x in [self.wav, self.val, self.unc])

    @property
    def ap_min(self):
        """
        Minimum aperture used to define the SEDs (in cm).
        """
        return self._ap_min

    @ap_min.setter
    def ap_min(self, value):
        if value is None or (np.isscalar(value) and np.isreal(value)):
            self._ap_min = value
        else:
            raise ValueError("ap_min should be a real scalar value")

    @property
    def ap_max(self):
        """
        Maximum aperture used to define the SEDs (in cm).
        """
        return self._ap_max

    @ap_max.setter
    def ap_max(self, value):
        if value is None or (np.isscalar(value) and np.isreal(value)):
            self._ap_max = value
        else:
            raise ValueError("ap_max should be a real scalar value")

    @property
    def d_min(self):
        """
        Minimum depth used to define the SEDs (in cm).
        """
        return self._d_min

    @d_min.setter
    def d_min(self, value):
        if value is None or (np.isscalar(value) and np.isreal(value)):
            self._d_min = value
        else:
            raise ValueError("d_min should be a real scalar value")

    @property
    def d_max(self):
        """
        Maximum depth used to define the SEDs (in cm).
        """
        return self._d_max

    @d_max.setter
    def d_max(self, value):
        if value is None or (np.isscalar(value) and np.isreal(value)):
            self._d_max = value
        else:
            raise ValueError("d_max should be a real scalar value")

    @property
    def distance(self):
        """
        Distance assumed for the image (in cm).
        """
        return self._distance

    @distance.setter
    def distance(self, value):
        if value is None or (np.isscalar(value) and np.isreal(value)):
            self._distance = value
        else:
            raise ValueError("distance should be a real scalar value")

    @property
    def inside_observer(self):
        """
        Whether the image was from an inside observer.
        """
        return self._inside_observer

    @inside_observer.setter
    def inside_observer(self, value):
        if value is None or type(value) is bool:
            self._inside_observer = value
        else:
            raise ValueError("inside_observer should be a boolean")
