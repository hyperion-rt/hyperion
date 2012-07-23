from __future__ import print_function, division

import numpy as np

from ._interpolate_core import interp1d_linear_scalar, \
                               interp1d_linear_array, \
                               interp1d_loglog_scalar, \
                               interp1d_loglog_array, \
                               interp1d_linlog_scalar, \
                               interp1d_linlog_array, \
                               interp1d_loglin_scalar, \
                               interp1d_loglin_array


class check_bounds(object):
    '''
    Decorator to add standard interpolation bounds checking.

    This decorator incurs an overhead of 30-50 percent on the runtime of
    the interpolation routine.
    '''

    def __init__(self, f):
        self.f = f

    def __call__(self, x, y, xval, bounds_error=True, fill_value=np.nan):
        if np.isscalar(xval):  # xval is a scalar
            if xval < x[0] or xval > x[-1]:  # the value is out of bounds
                if bounds_error:
                    raise Exception("x value is out of interpolation bounds")
                else:
                    return fill_value
            else:  # the value is in the bounds
                return self.f(x, y, xval)
        else:  # xval is an array
            inside = (xval >= x[0]) & (xval <= x[-1])
            outside = ~inside
            if np.any(outside):  # some values are out of bounds
                if bounds_error:
                    raise Exception("x values are out of interpolation bounds")
                else:
                    if np.any(inside):
                        yval = np.zeros(xval.shape)
                        yval[inside] = self.f(x, y, xval[inside])
                        yval[outside] = fill_value
                        return yval
                    else:
                        return np.repeat(fill_value, xval.shape)
            else:  # all values are in the bounds
                return self.f(x, y, xval)


@check_bounds
def interp1d_fast(x, y, xval):
    '''On-the-fly linear interpolator'''
    if len(x) != len(y):
        raise Exception("x and y should have the same length")
    if x.dtype != float or y.dtype != float:
        x, y = x.astype(float), y.astype(float)
    if np.isscalar(xval):
        return interp1d_linear_scalar(x, y, np.float(xval))
    else:
        if xval.ndim > 1:
            return interp1d_linear_array(x, y, xval.flatten()).reshape(xval.shape)
        else:
            return interp1d_linear_array(x, y, xval)


@check_bounds
def interp1d_fast_loglog(x, y, xval):
    '''On-the-fly log interpolator'''
    if len(x) != len(y):
        raise Exception("x and y should have the same length")
    if x.dtype != float or y.dtype != float:
        x, y = x.astype(float), y.astype(float)
    if np.isscalar(xval):
        return interp1d_loglog_scalar(x, y, np.float(xval))
    else:
        if xval.ndim > 1:
            return interp1d_loglog_array(x, y, xval.flatten()).reshape(xval.shape)
        else:
            return interp1d_loglog_array(x, y, xval)


@check_bounds
def interp1d_fast_linlog(x, y, xval):
    '''On-the-fly linear-log interpolator'''
    if len(x) != len(y):
        raise Exception("x and y should have the same length")
    if x.dtype != float or y.dtype != float:
        x, y = x.astype(float), y.astype(float)
    if np.isscalar(xval):
        return interp1d_linlog_scalar(x, y, np.float(xval))
    else:
        if xval.ndim > 1:
            return interp1d_linlog_array(x, y, xval.flatten()).reshape(xval.shape)
        else:
            return interp1d_linlog_array(x, y, xval)


@check_bounds
def interp1d_fast_loglin(x, y, xval):
    '''On-the-fly log-linear interpolator'''
    if len(x) != len(y):
        raise Exception("x and y should have the same length")
    if x.dtype != float or y.dtype != float:
        x, y = x.astype(float), y.astype(float)
    if np.isscalar(xval):
        return interp1d_loglin_scalar(x, y, np.float(xval))
    else:
        if xval.ndim > 1:
            return interp1d_loglin_array(x, y, xval.flatten()).reshape(xval.shape)
        else:
            return interp1d_loglin_array(x, y, xval)
