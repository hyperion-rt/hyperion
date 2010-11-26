import warnings

import numpy as np

from hyperion.util.constants import pi
from hyperion.util.functions import FreezableClass
from hyperion.densities.bipolar_cavity import BipolarCavity
from hyperion.util.convenience import OptThinRadius


class PowerLawEnvelope(FreezableClass):

    def __init__(self, mass=None, rmin=None, rmax=None, power=None):
        '''
        Initialize a power-law envelope instance. The required parameters are:

            mass: envelope mass (g)
            rmin: inner radius (cm)
            rmin: output radius (cm)
            power: power-law exponent

        '''

        # Basic envelope parameters
        self.mass = mass
        self.rmin = rmin
        self.rmax = rmax
        self.power = power

        # Cavity
        self.cavity = None

        # Dust
        self.dust = None

        self._freeze()

    def _check_all_set(self):

        if self.mass is None:
            raise Exception("mass is not set")
        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")
        if self.power is None:
            raise Exception("power is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner envelope radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Outer envelope radius needs to be computed first")

    def exists(self):
        return self.mass > 0.

    def density(self, grid):
        '''
        Find the density of a power-law envelope

        Input is the position of the center of the grid cells in spherical
        polar coordinates (r, theta, phi), the volume of the grid cells, and a
        parameter dictionary
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring power-law envelope, since rmax < rmin")
            return np.zeros(grid.shape)

        alpha = 3. + self.power

        rho_0 = self.mass * alpha / \
            (4. * pi * (self.rmax ** alpha - self.rmin ** alpha))

        rho = rho_0 * grid.gr ** self.power

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        norm = self.mass / np.sum(rho * grid.volumes)

        print "Normalization factor for envelope mass (should ideally be 1): %5.2f" % norm

        rho = rho * norm

        if self.cavity is not None:
            mask = self.cavity.mask(grid)
            rho[~mask] = 0.

        return rho

    def midplane_cumulative_density(self, r):
        '''
        Find the cumulative column density as a function of radius from the
        star in the midplane of a power-law envelope.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring power-law envelope, since rmax < rmin")
            return np.zeros(r.shape)

        p = 1 - self.power

        if p == 0.:
            rho = self.rho_0 * np.log(r / self.rmin)
        else:
            rho = self.rho_0 * (r**p - self.rmin**p) / p

        return rho

    def add_bipolar_cavity(self):
        if self.cavity is not None:
            raise Exception("Envelope already has a bipolar cavity")
        else:
            self.cavity = BipolarCavity()
            return self.cavity
