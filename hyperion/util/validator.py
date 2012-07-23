import numpy as np


def validate_scalar(name, value, domain=None, extra=''):

    if not np.isscalar(value):
        raise ValueError("{0:s} should be a scalar value{1:s}".format(name, extra))
    if not np.isreal(value):
        raise ValueError("{0:s} should be a numerical value{1:s}".format(name, extra))

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
