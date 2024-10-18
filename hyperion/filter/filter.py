from __future__ import print_function, division

import numpy as np

import six
from astropy import units as u

from ..util.integrate import integrate
from ..util.validator import validate_scalar, validate_array


class Filter(object):
    """
    A 'filter', in the general sense of a spectral transmission curve.

    The filter should be specified as an absolute transmission between 0 and
    100% as a function of frequency. This will then be used to compute a total
    energy in ergs/s or ergs/s/cm^2. It is left to the user to decide how to
    convert this to a monochromatic frequency.

    Parameters
    ----------
    name : str
        The name of the filter.
    spectral_coord : :class:`~astropy.units.quantity.Quantity`
        The spectral coordinates (e.g. wavelength, frequency, photon energy) at
        which the transmissions are defined.
    transmission : `numpy.ndarray`
        The spectral transmission as a function of spectral coordinate.
    """

    def __init__(self, name=None, spectral_coord=None, transmission=None):
        self.name = name
        self.spectral_coord = spectral_coord
        self.transmission = transmission
        self._alpha = None
        self._beta = None
        self.central_spectral_coord = None

    @property
    def name(self):
        """
        The name of the filter.
        """
        return self._name

    @name.setter
    def name(self, value):
        if value is None or isinstance(value, six.string_types):
            self._name = value
        else:
            raise TypeError("name should be given as a string")

    @property
    def spectral_coord(self):
        """
        The spectral coordinates (e.g. wavelength, frequency, photon energy) at
        which the filter is defined.
        """
        return self._spectral_coord

    @spectral_coord.setter
    def spectral_coord(self, value):
        if value is None:
            self._spectral_coord = None
        else:
            self._spectral_coord = validate_array('spectral_coord', value, domain='strictly-positive', ndim=1,
                                                  physical_type=('frequency', 'length', 'energy'))

    @property
    def transmission(self):
        """
        The filter transmission.
        """
        return self._transmission

    @transmission.setter
    def transmission(self, value):
        if value is None:
            self._transmission = None
        else:
            self._transmission = validate_array('r', value, domain='positive', ndim=1,
                                                shape=None if self.spectral_coord is None else (len(self.spectral_coord),),
                                                physical_type=('dimensionless'))

    def check_all_set(self):
        for attr in ['spectral_coord', 'transmission', 'name', 'alpha',
                     'detector_type', 'central_spectral_coord']:
            if getattr(self, attr) is None:
                raise ValueError("{0} has not been set".format(attr))

    def to_hdf5_group(self, group, name):

        self.check_all_set()

        # Get spectral coordinate in Hz and transmision in fractional terms
        nu = self.spectral_coord.to(u.Hz, equivalencies=u.spectral()).value
        tr = self.transmission.to(u.one).value

        # Sort in order of increasing Hz
        order = np.argsort(nu)
        nu = nu[order]
        tr = tr[order]

        # Get other parameters for the normalization
        nu0 = self.central_spectral_coord.to(u.Hz, equivalencies=u.spectral()).value
        alpha = self.alpha
        beta = self._beta

        # Here we normalize the filter before passing it to Hyperion
        tr_norm = (tr / nu ** (1 + beta)
                   / nu0 ** alpha
                   / integrate(nu, tr / nu ** (1. + alpha + beta)))

        # Now multiply by nu so that Hyperion returns nu * Fnu
        tr_norm *= nu

        dset = group.create_dataset(name, data=np.array(list(zip(nu, tr, tr_norm)),
                                                        dtype=[('nu', float),
                                                               ('tr', float),
                                                               ('tn', float)]))

        dset.attrs['name'] = np.bytes_(self.name)

        dset.attrs['alpha'] = self.alpha
        dset.attrs['beta'] = self._beta
        dset.attrs['nu0'] = self.central_spectral_coord.to(u.Hz, equivalencies=u.spectral()).value

    @classmethod
    def from_hdf5_group(cls, group, name):

        self = cls()
        self.spectral_coord = group[name]['nu'] * u.Hz
        self.transmission = group[name]['tr'] * u.one
        self.name = group[name].attrs['name'].decode('utf-8')

        self.alpha = group[name].attrs['alpha']
        self._beta = group[name].attrs['beta']
        self.central_spectral_coords = group[name].attrs['nu0'] * u.Hz

        return self

    @property
    def detector_type(self):
        return "energy" if self._beta == -1 else "photons"

    @detector_type.setter
    def detector_type(self, value):
        if value == 'energy':
            self._beta = -1
        elif value == 'photons':
            self._beta = 0
        else:
            raise ValueError("detector_type should be one of energy/photons")

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @property
    def central_spectral_coord(self):
        """
        The central spectral coordinate (e.g. wavelength, frequency, photon energy) at
        which the monochromatic flux should be measured.
        """
        return self._central_spectral_coord

    @central_spectral_coord.setter
    def central_spectral_coord(self, value):
        if value is None:
            self._central_spectral_coord = None
        else:
            self._central_spectral_coord = validate_scalar('spectral_coord', value, domain='strictly-positive',
                                                           physical_type=('frequency', 'length', 'energy'))
