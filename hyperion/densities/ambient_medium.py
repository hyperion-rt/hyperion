import numpy as np

from hyperion.util.functions import FreezableClass
from hyperion.util.convenience import OptThinRadius


class AmbientMedium(FreezableClass):

    def __init__(self, density=None, temperature=2.725, rmin=None):
        "Initialize an ambient medium instance"

        # Basic ambient medium parameters
        self.density = density
        self.temperature = temperature
        self.rmin = rmin
        self.rmax = rmax

        # Dust
        self.dust = None

        self._freeze()

    def check_all_set(self):

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

        self.check_all_set()

        rho = np.repeat(self.density, grid.gr.shape)

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        return rho
