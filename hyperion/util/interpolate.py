import numpy as np


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
                    return np.nan
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
    ipos = np.searchsorted(x, xval)
    return (xval - x[ipos - 1]) \
         / (x[ipos] - x[ipos - 1]) \
         * (y[ipos] - y[ipos - 1]) \
         + y[ipos - 1]


@check_bounds
def interp1d_fast_loglog(x, y, xval):
    '''On-the-fly log interpolator'''
    ipos = np.searchsorted(x, xval)
    yval = 10.**((np.log10(xval) - np.log10(x[ipos - 1])) \
              / (np.log10(x[ipos]) - np.log10(x[ipos - 1])) \
              * (np.log10(y[ipos]) - np.log10(y[ipos - 1])) \
              + np.log10(y[ipos - 1]))
    if np.isscalar(xval):
        if y[ipos - 1] == 0. or y[ipos] == 0.:
            yval = 0.
    else:
        yval[(y[ipos - 1] == 0.) | (y[ipos] == 0.)] = 0.
    return yval


@check_bounds
def interp1d_fast_linlog(x, y, xval):
    '''On-the-fly linear-log interpolator'''
    ipos = np.searchsorted(x, xval)
    yval = 10.**((xval - x[ipos - 1]) \
              / (x[ipos] - x[ipos - 1]) \
              * (np.log10(y[ipos]) - np.log10(y[ipos - 1])) \
              + np.log10(y[ipos - 1]))
    if np.isscalar(xval):
        if y[ipos - 1] == 0. or y[ipos] == 0.:
            yval = 0.
    else:
        yval[(y[ipos - 1] == 0.) | (y[ipos] == 0.)] = 0.
    return yval


@check_bounds
def interp1d_fast_loglin(x, y, xval):
    '''On-the-fly log-linear interpolator'''
    ipos = np.searchsorted(x, xval)
    return ((np.log10(xval) - np.log10(x[ipos-1])) \
           / (np.log10(x[ipos]) - np.log10(x[ipos - 1])) \
           * (y[ipos] - y[ipos - 1]) \
           + y[ipos - 1])

# Legacy

from scipy.interpolate import interp1d


def monotonic(array):
    for i in range(len(array)-1):
        if array[i+1] < array[i]:
            return False
    return True


class interp1d_log10(object):

    def __init__(self, x, y, bounds_error=True, fill_value=np.nan, **kwargs):

        if not monotonic(x):
            raise Exception("Array is not monotonically increasing")

        self.bounds_error = bounds_error
        self.fill_value = fill_value
        self.xmin = np.min(x)
        self.xmax = np.max(x)
        self.ymin = y[x==self.xmin]
        self.ymax = y[x==self.xmax]
        self.interp = interp1d(np.log10(x), np.log10(y), **kwargs)
        self.interp_lin = interp1d(x, y, **kwargs)

    def __call__(self, x_new):

        if type(x_new) == np.ndarray:

            if self.bounds_error:
                if np.any((x_new < self.xmin) | (x_new > self.xmax)):
                    string = "\n\n"
                    for x in x_new.ravel():
                        if x < self.xmin or x > self.xmax:
                            string += "    xmin = %g\n" % self.xmin
                            string += "    x    = %g\n" % x
                            string += "    xmax = %g\n" % self.xmax
                        raise Exception("Interpolation out of bounds:" + string)

            y_new = np.ones(x_new.shape) * self.fill_value

            y_new[x_new == self.xmin] = self.ymin
            y_new[x_new == self.xmax] = self.ymax

            mask = (x_new >= self.xmin) & (x_new <= self.xmax)
            y_new[mask] = 10.**self.interp(np.log10(x_new[mask]))
            y_new[mask & np.isnan(y_new)] = self.interp_lin(x_new[mask & np.isnan(y_new)])

        else:

            if x_new < self.xmin or x_new > self.xmax:
                if self.bounds_error:
                    raise Exception("Interpolation out of bounds")
                else:
                    y_new = self.fill_value
            elif x_new == self.xmin:
                y_new = self.ymin
            elif x_new == self.xmax:
                y_new = self.ymax
            else:
                y_new = 10.**self.interp(np.log10(x_new))
                if np.isnan(y_new):
                    y_new = self.interp_lin(x_new)

        return y_new
