import warnings

import numpy as np

from hyperion.util.constants import pi, G
from hyperion.util.functions import FreezableClass
from hyperion.util.convenience import OptThinRadius
from hyperion.util.integrate import integrate_powerlaw
from hyperion.dust import SphericalDust


class FlaredDisk(FreezableClass):

    def __init__(self):
        '''
        Initialize a flared disk instance.

        The available attributes are:

            mass: mass (g)
            rmin: inner radius (cm)
            rmax: outer radius (cm)
            p: surface density exponent
            beta: flaring exponent
            h_0: scaleheight at r_0 (cm)
            r_0: radius at which scaleheight is defined (cm)

        '''

        # Basic disk parameters
        self.mass = 0.
        self.rmin = None
        self.rmax = None
        self.p = -1.
        self.beta = 1.25
        self.h_0 = None
        self.r_0 = None

        # Fine control over density distribution
        self.cylindrical_inner_rim = True
        self.cylindrical_outer_rim = True

        # Dust
        self.dust = None

        self._freeze()

    def __str__(self):
        string = "= Flared disk =\n"
        string += " - M_disk: %.3e\n" % self.mass
        string += " - R_min: %.3e\n" % self.rmin
        string += " - R_min: %.3e\n" % self.rmax
        string += " - p: %.3f\n" % self.p
        string += " - beta: %.3f\n" % self.beta
        string += " - h_0: %.3e\n" % self.h_0
        string += " - r_0: %.3e\n" % self.r_0
        return string

    def _check_all_set(self):

        if self.mass is None:
            raise Exception("mass is not set")
        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")
        if self.p is None:
            raise Exception("p is not set")
        if self.beta is None:
            raise Exception("beta is not set")
        if self.h_0 is None:
            raise Exception("h_0 is not set")
        if self.r_0 is None:
            raise Exception("r_0 is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner disk radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Outer disk radius needs to be computed first")

    def rho_0(self):
        '''
        Returns the density factor rho0
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring disk, since rmax < rmin")
            return 0.

        int1 = integrate_powerlaw(self.rmin, self.rmax, 1.0 + self.p)
        int1 *= self.r_0 ** -self.p

        integral = (2. * pi) ** 1.5 * self.h_0 * int1

        return self.mass / integral

    def density(self, grid):
        '''
        Return a density grid for spherical polar coordinates

        Input is the position of the center of the grid cells in spherical
        polar coordinates (r, theta, phi), and the volume of the grid cells.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(grid.shape)

        if self.mass == 0:
            return np.zeros(grid.shape)

        # Find disk scaleheight at each cylindrical radius
        h = self.h_0 * (grid.gw / self.r_0) ** self.beta

        # Find disk density at all positions
        rho = (self.r_0 / grid.gw) ** (self.beta - self.p) \
            * np.exp(-0.5 * (grid.gz / h) ** 2)

        # Truncate below rmin and above rmax
        if self.cylindrical_inner_rim:
            rho[grid.gw < self.rmin] = 0.
        else:
            rho[grid.gr < self.rmin] = 0.

        if self.cylindrical_outer_rim:
            rho[grid.gw > self.rmax] = 0.
        else:
            rho[grid.gr > self.rmax] = 0.

        # Find density factor
        rho *= self.rho_0()

        norm = self.mass / np.sum(rho * grid.volumes)

        print "Normalization factor for disk mass: %5.2f" % norm

        # Normalize to total disk mass
        rho = rho * norm

        return rho

    def midplane_cumulative_density(self, r):
        '''
        Find the cumulative column density as a function of radius from the
        star in the midplane of a standard flared disk.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(r.shape)

        int1 = integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), self.p - self.beta)
        int1 *= self.r_0 ** (self.beta - self.p)

        return self.rho_0() * int1

    def vertical_profile(self, r, theta):

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(theta.shape)

        # Convert coordinates to cylindrical polars
        z = r * np.cos(theta)
        w = r * np.sin(theta)

        # Find disk scaleheight at each cylindrical radius
        h = self.h_0 * (w / self.r_0) ** self.beta

        # Find disk density at all positions
        rho = (self.r_0 / w) ** (self.beta - self.p) \
            * np.exp(-0.5 * (z / h) ** 2)

        rho *= self.rho_0()

        # What about normalization

        return rho

    def vertical_cumulative_density(self, r, theta):

        density = self.vertical_profile(r, theta)

        d = r * np.radians(theta)

        tau = density * d

        tau[0] = 0.

        return tau

    def __setattr__(self, attribute, value):
        if attribute == 'dust' and value is not None and type(value) is str:
            FreezableClass.__setattr__(self, 'dust', SphericalDust(value))
        else:
            FreezableClass.__setattr__(self, attribute, value)
