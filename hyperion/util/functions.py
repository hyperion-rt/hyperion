import string
import random
import os
import shutil
import sys

import numpy as np

from hyperion.util.constants import h, c, k


class FreezableClass(object):

    _frozen = False
    _final = False

    def _freeze(self):
        object.__setattr__(self, '_frozen', True)

    def _finalize(self):
        object.__setattr__(self, '_final', True)

    def isfrozen(self):
        return self._frozen

    def isfinal(self):
        return self._final

    def __setattr__(self, key, value):
        if self._final:
            raise Exception("Attribute %s can no longer be changed" % key)
        if self._frozen and not hasattr(self, key):
            raise AttributeError("Attribute %s does not exist" % key)
        object.__setattr__(self, key, value)


def planck_nu_range(tmin, tmax=None):

    # Set constant for Wien displacement law
    alpha = 2.821439

    # Find peak frequencies for lower and upper temperature
    nu_peak_min = alpha / h * k * tmin
    if tmax is None:
        nu_peak_max = alpha / h * k * tmin
    else:
        nu_peak_max = alpha / h * k * tmax

    # Find frequency range if we want to go below 1 thousandth of the Planck
    # function peak

    nu_min = np.log10(nu_peak_min / 100.)
    nu_max = np.log10(nu_peak_max * 10.)

    # If we use at least 100 frequencies per order of magnitude of frequency
    # then we will achieve differences around 2% flux-to-flux

    n_nu = int((nu_max - nu_min) * 100.)

    return np.logspace(nu_min, nu_max, n_nu)


def nu_common(nu1, nu2):

    # Combine frequency lists/arrays and sort
    nu_common = np.hstack([nu1, nu2])
    nu_common.sort()

    # Remove unique elements (can't just use np.unique because also want
    # to remove very close values)
    keep = (nu_common[1:] - nu_common[:-1]) / nu_common[:-1] > 1.e-10
    keep = np.hstack([keep, True])
    nu_common = nu_common[keep]

    # Return common frequency range
    return nu_common


class extrap1d_log10(object):

    def __init__(self, x, y):

        self.x = np.log10(x)
        self.y = np.log10(y)

    def __call__(self, x):

        xval = np.log10(x)

        if self.x[-1] > self.x[0]:
            top = xval > self.x[-1]
            bot = xval < self.x[0]
        else:
            top = xval < self.x[-1]
            bot = xval > self.x[0]

        if type(xval) == np.ndarray:

            yval = np.zeros(xval.shape)

            yval[top] = 10. ** (self.y[-1] + (xval[top] - self.x[-1]) * (self.y[-1] - self.y[-2]) / (self.x[-1] - self.x[-2]))
            yval[bot] = 10. ** (self.y[0] + (xval[bot] - self.x[0]) * (self.y[0] - self.y[1]) / (self.x[0] - self.x[1]))

        else:

            if top:
                yval = 10. ** (self.y[-1] + (xval - self.x[-1]) * (self.y[-1] - self.y[-2]) / (self.x[-1] - self.x[-2]))
            elif bot:
                yval = 10. ** (self.y[0] + (xval - self.x[0]) * (self.y[0] - self.y[1]) / (self.x[0] - self.x[1]))
            else:
                raise Exception("xval should lie outside x array")

        return yval


def B_nu(nu, T):
    return 2. * h * nu ** 3. / c ** 2. / (np.exp(h * nu / k / T) - 1.)


def filename2fits(filename):

    ext = os.path.splitext(filename)[1]
    if ext in ['.par', '.dat', '.txt', '.ascii']:
        return filename.replace(ext, '.fits')
    else:
        raise Exception("Unknown extension: %s" % ext)


def filename2hdf5(filename):

    ext = os.path.splitext(filename)[1]
    if ext in ['.par', '.dat', '.txt', '.ascii']:
        return filename.replace(ext, '.hdf5')
    else:
        raise Exception("Unknown extension: %s" % ext)


def random_id(length=32):
    chars = string.letters + string.digits
    s = ""
    for i in range(length):
        s += random.choice(chars)
    return s


def create_dir(dir_name):
    delete_dir(dir_name)
    os.mkdir(dir_name)


def delete_dir(dir_name):

    if os.path.exists(dir_name):
        reply = raw_input("Delete directory %s? [y/[n]] " % dir_name)
        if reply == 'y':
            shutil.rmtree(dir_name)
        else:
            print "Aborting..."
            sys.exit()


def delete_file(file_name):

    if os.path.exists(file_name):
        reply = raw_input("Delete file %s? [y/[n]] " % file_name)
        if reply == 'y':
            os.remove(file_name)
        else:
            print "Aborting..."
            sys.exit()
