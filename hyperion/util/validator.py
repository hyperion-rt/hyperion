import numpy as np

from astropy import units as u
from astropy.extern import six


def validate_physical_type(name, value, physical_type):
    if physical_type is not None:
        if not isinstance(value, u.Quantity):
            raise TypeError("{0} should be given as a Quantity object".format(name))
        if isinstance(physical_type, six.string_types):
            if value.unit.physical_type != physical_type:
                raise TypeError("{0} should be given in units of {1}".format(name, physical_type))
        else:
            if not value.unit.physical_type in physical_type:
                raise TypeError("{0} should be given in units of {1}".format(name, ', '.join(physical_type)))


def validate_array(name, value, domain=None, ndim=1, shape=None, physical_type=None):

    if physical_type is None:
        if type(value) in [list, tuple]:
            value = np.array(value)
    else:
        validate_physical_type(name, value, physical_type)

    # Check the value is an array with the right number of dimensions
    if not isinstance(value, np.ndarray) or value.ndim != ndim:
        if ndim == 1:
            raise TypeError("{0} should be a 1-d sequence".format(name))
        else:
            raise TypeError("{0} should be a {1:d}-d array".format(name, ndim))

    # Check that the shape matches that expected
    if shape is not None and value.shape != shape:
        if ndim == 1:
            raise ValueError("{0} has incorrect length (expected {1} but found {2})".format(name, shape[0], value.shape[0]))
        else:
            raise ValueError("{0} has incorrect shape (expected {1} but found {2})".format(name, shape, value.shape))

    return value


def validate_scalar(name, value, domain=None, extra='', physical_type=None):

    if physical_type is None:

        if not np.isscalar(value):
            raise ValueError("{0:s} should be a scalar value{1:s}".format(name, extra))

        if not np.isreal(value):
            raise ValueError("{0:s} should be a numerical value{1:s}".format(name, extra))

    else:

        validate_physical_type(name, value, physical_type)

    if domain == 'positive':
        if value < 0.:
            raise ValueError("{0:s} should be positive".format(name))
    elif domain == 'strictly-positive':
        if value <= 0.:
            raise ValueError("{0:s} should be strictly positive".format(name))
    elif domain == 'negative':
        if value > 0.:
            raise ValueError("{0:s} should be negative".format(name))
    elif domain == 'strictly-negative':
        if value >= 0.:
            raise ValueError("{0:s} should be strictly negative".format(name))
    elif type(domain) in [tuple, list] and len(domain) == 2:
        if value < domain[0] or value > domain[-1]:
            raise ValueError("{0:s} should be in the range [{1:g}:{2:g}]".format(name, domain[0], domain[-1]))

    return value
