import numpy as np


def interp1d_fast(x, y, xval):
    '''On-the-fly linear interpolator'''
    ipos = np.searchsorted(x, xval)
    return (xval - x[ipos - 1]) \
         / (x[ipos] - x[ipos - 1]) \
         * (y[ipos] - y[ipos - 1]) \
         + y[ipos - 1]


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
