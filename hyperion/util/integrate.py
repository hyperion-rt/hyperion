import numpy as np
from interpolate import interp1d_fast, interp1d_fast_loglog, \
                        interp1d_fast_linlog, interp1d_fast_loglin


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
        ymin = interp1d_fast(x[i1-1:i1+1], y[i1-1:i1+1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast(x[i2-1:i2+1], y[i2-1:i2+1], xmax)

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
        ymin = interp1d_fast_loglin(x[i1-1:i1+1], y[i1-1:i1+1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_loglin(x[i2-1:i2+1], y[i2-1:i2+1], xmax)

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
        ymin = interp1d_fast_linlog(x[i1-1:i1+1], y[i1-1:i1+1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_linlog(x[i2-1:i2+1], y[i2-1:i2+1], xmax)

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
        ymin = interp1d_fast_loglog(x[i1-1:i1+1], y[i1-1:i1+1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast_loglog(x[i2-1:i2+1], y[i2-1:i2+1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate_loglog(x, y)


def integrate(x, y):
    '''
    Perform trapezium integration of a set of points (x,y). The interpolation
    between the points is done in linear space, so this is designed for
    functions that are piecewise linear in linear space.
    '''

    # Fix NaN values
    y[np.isnan(y)] = 0.

    # Find the integral of all the chunks
    integrals = 0.5 * (x[1:] - x[:-1]) * (y[1:] + y[:-1])

    # Sum them all up
    integral = np.sum(integrals)

    # Check if the integral is NaN or infinity
    if np.isnan(integral) or np.isinf(integral):
        raise Exception("Integral is NaN or Inf")

    return integral


def integrate_loglin(x, y):
    '''
    Perform trapezium integration of a set of points (x,y). The interpolation
    between the points is done in log-linear space, so this is designed for
    functions that are piecewise linear in log-linear space.
    '''

    # Fix NaN values
    y[np.isnan(y)] = 0.

    # Compute the power of the power-laws connecting the points
    a = (y[:-1] - y[1:]) / np.log10(x[:-1] / x[1:])
    b = y[:-1] - a * np.log10(x[:-1])

    # Find the integral of all the chunks
    integrals = a * (x[1:] * np.log10(x[1:]) - x[:-1] * np.log10(x[:-1])) \
              + (b - a / np.log(10.)) * (x[1:] - x[:-1])

    reset = x[1:] == x[:-1]
    if np.any(reset):
        integrals[reset] = 0.

    # Sum them all up
    integral = np.sum(integrals)

    # Check if the integral is NaN or infinity
    if np.isnan(integral) or np.isinf(integral):
        raise Exception("Integral is NaN or Inf")

    return integral


def integrate_linlog(x, y):
    '''
    Perform trapezium integration of a set of points (x,y). The interpolation
    between the points is done in linear-log space, so this is designed for
    functions that are piecewise linear in linear-log space.
    '''

    # Fix NaN values
    y[np.isnan(y)] = 0.

    # Find which bins to ignore
    keep = (y[:-1] > 0.) & (y[1:] > 0.)

    # Find the integral of all the chunks
    integrals = (y[1:] - y[:-1]) * (x[1:] - x[:-1]) \
              / np.log(10.) / np.log10(y[1:] / y[:-1])

    reset = x[1:] == x[:-1]
    if np.any(reset):
        integrals[reset] = 0.

    reset = y[1:] == y[:-1]
    if np.any(reset):
        integrals[reset] = y[:-1] * (x[1:] - x[:-1])

    # Sum them all up
    integral = np.sum(integrals[keep])

    # Check if the integral is NaN or infinity
    if np.isnan(integral) or np.isinf(integral):
        raise Exception("Integral is NaN or Inf")

    return integral


def integrate_loglog(x, y):
    '''
    Perform trapezium integration of a set of points (x,y). The interpolation
    between the points is done in log-log space, so this is designed for
    functions that are piecewise linear in log-log space.
    '''

    # Fix NaN values
    y[np.isnan(y)] = 0.

    # Find which bins to ignore
    keep = (y[:-1] > 0.) & (y[1:] > 0.)

    # Compute the power of the power-laws connecting the points
    b = np.log10(y[:-1] / y[1:]) / np.log10(x[:-1] / x[1:])

    # Find the integral of all the chunks
    integrals = y[:-1] * (x[1:] * np.power(x[1:] / x[:-1], b) - x[:-1]) \
              / (b + 1)

    # Address special case
    reset = np.abs(b + 1.) < 1e-10
    if np.any(reset):
        integrals[reset] = x[:-1][reset] * y[:-1][reset] \
                         * np.log(x[1:][reset] / x[:-1][reset])

    # Sum them all up
    integral = np.sum(integrals[keep])

    # Check if the integral is NaN or infinity
    if np.isnan(integral) or np.isinf(integral):
        raise Exception("Integral is NaN or Inf")

    return integral


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
