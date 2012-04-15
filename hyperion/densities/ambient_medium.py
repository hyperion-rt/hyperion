from __future__ import print_function, division

import numpy as np

from hyperion.util.functions import FreezableClass
from hyperion.util.convenience import OptThinRadius
from hyperion.dust import SphericalDust


class AmbientMedium(FreezableClass):

    def __init__(self, rho=None, temperature=2.725, rmin=None, rmax=None):
        "Initialize an ambient medium instance"

        # Basic ambient medium parameters
        self.rho = rho
        self.temperature = temperature
        self.rmin = rmin
        self.rmax = rmax

        # Dust
        self.dust = None

        self._freeze()

    def _check_all_set(self):

        if self.density is None:
            raise Exception("density is not set")
        if self.temperature is None:
            raise Exception("temperature is not set")
        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner ambient medium radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Inner ambient medium radius needs to be computed first")

    def density(self, grid):
        '''
        Find the density of the ambient medium

        Input is the position of the center of the grid cells in spherical
        polar coordinates (r, theta, phi), the volume of the grid cells, and a
        parameter dictionary
        '''

        self._check_all_set()

        rho = np.ones(grid.gr.shape) * self.rho

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        return rho

    def __setattr__(self, attribute, value):
        if attribute == 'dust' and value is not None and type(value) is str:
            FreezableClass.__setattr__(self, 'dust', SphericalDust(value))
        else:
            FreezableClass.__setattr__(self, attribute, value)
