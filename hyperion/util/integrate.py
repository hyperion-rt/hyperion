from __future__ import print_function, division

import numpy as np

from .interpolate import interp1d_fast, \
                         interp1d_fast_loglin, \
                         interp1d_fast_linlog, \
                         interp1d_fast_loglog

from ._integrate_core import _integrate, \
                             _integrate_loglin, \
                             _integrate_linlog, \
                             _integrate_loglog


def integrate_subset(x, y, xmin, xmax):
    '''
    Perform trapezium integration of a set of points (x,y) between bounds xmin
    and xmax. The interpolation between the points is done in linear space, so
    this is designed for functions that are piecewise linear in linear space.
    '''

    # Swap arrays if necessary
    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    # Swap limits if necessary
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    elif xmin == xmax:
        return 0.

    # Find the subset of points to use and the value of the function at the
    # end-points of the integration

    if xmin == x[0]:
        i1 = 1
        ymin = y[0]
    else:
        i1 = np.searchsorted(x, xmin)
        if xmin == x[i1]:
            i1 += 1
        ymin = interp1d_fast(x[i1 - 1:i1 + 1], y[i1 - 1:i1 + 1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast(x[i2 - 1:i2 + 1], y[i2 - 1:i2 + 1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate(x, y)


def integrate_loglin_subset(x, y, xmin, xmax):
    '''
    Perform trapezium integration of a set of points (x,y) between bounds xmin
    and xmax. The interpolation between the points is done in log-linear
    space, so this is designed for functions that are piecewise linear in
    log-linear space.
    '''

    # Swap arrays if necessary
    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    # Swap limits if necessary
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    elif xmin == xmax:
        return 0.

    # Find the subset of points to use and the value of the function at the
    # end-points of the integration

    if xmin == x[0]:
        i1 = 1
        ymin = y[0]
    else:
        i1 = np.searchsorted(x, xmin)
        if xmin == x[i1]:
            i1 += 1
        ymin = interp1d_fast_loglin(x[i1 - 1:i1 + 1], y[i1 - 1:i1 + 1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_loglin(x[i2 - 1:i2 + 1], y[i2 - 1:i2 + 1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate_loglin(x, y)


def integrate_linlog_subset(x, y, xmin, xmax):
    '''
    Perform trapezium integration of a set of points (x,y) between bounds xmin
    and xmax. The interpolation between the points is done in linear-log
    space, so this is designed for functions that are piecewise linear in
    linear-log space.
    '''

    # Swap arrays if necessary
    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    # Swap limits if necessary
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    elif xmin == xmax:
        return 0.

    # Find the subset of points to use and the value of the function at the
    # end-points of the integration

    if xmin == x[0]:
        i1 = 1
        ymin = y[0]
    else:
        i1 = np.searchsorted(x, xmin)
        if xmin == x[i1]:
            i1 += 1
        ymin = interp1d_fast_linlog(x[i1 - 1:i1 + 1], y[i1 - 1:i1 + 1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_linlog(x[i2 - 1:i2 + 1], y[i2 - 1:i2 + 1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate_linlog(x, y)


def integrate_loglog_subset(x, y, xmin, xmax):
    '''
    Perform trapezium integration of a set of points (x,y) between bounds xmin
    and xmax. The interpolation between the points is done in log-log space, so
    this is designed for functions that are piecewise linear in log-log space.
    '''

    # Swap arrays if necessary
    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    # Swap limits if necessary
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    elif xmin == xmax:
        return 0.

    # Find the subset of points to use and the value of the function at the
    # end-points of the integration

    if xmin == x[0]:
        i1 = 1
        ymin = y[0]
    else:
        i1 = np.searchsorted(x, xmin)
        if xmin == x[i1]:
            i1 += 1
        ymin = interp1d_fast_loglog(x[i1 - 1:i1 + 1], y[i1 - 1:i1 + 1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_loglog(x[i2 - 1:i2 + 1], y[i2 - 1:i2 + 1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate_loglog(x, y)


def integrate(x, y):
    if x.dtype == float and y.dtype == float:
        return _integrate(x, y)
    else:
        return _integrate(x.astype(float), y.astype(float))


def integrate_loglin(x, y):
    if x.dtype == float and y.dtype == float:
        return _integrate_loglin(x, y)
    else:
        return _integrate_loglin(x.astype(float), y.astype(float))


def integrate_linlog(x, y):
    if x.dtype == float and y.dtype == float:
        return _integrate_linlog(x, y)
    else:
        return _integrate_linlog(x.astype(float), y.astype(float))


def integrate_loglog(x, y):
    if x.dtype == float and y.dtype == float:
        return _integrate_loglog(x, y)
    else:
        return _integrate_loglog(x.astype(float), y.astype(float))


def integrate_powerlaw(xmin, xmax, power):
    '''
    Find the integral of:

         xmax
        /
        | x^power dx
        /
     xmin
    '''
    if power == -1.:
        return np.log(xmax / xmin)
    else:
        return (xmax ** (power + 1.) - xmin ** (power + 1.)) / (power + 1.)
