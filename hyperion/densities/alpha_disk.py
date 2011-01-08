import warnings

import numpy as np

from hyperion.util.constants import pi, G
from hyperion.util.functions import FreezableClass
from hyperion.util.convenience import OptThinRadius
from hyperion.util.integrate import integrate_powerlaw
from hyperion.dust import SphericalDust
from hyperion.densities.flared_disk import FlaredDisk


class AlphaDisk(FlaredDisk):

    def __init__(self, mdot=None, lvisc=None, star=None, **kwargs):
        '''
        Initialize a flared disk instance. The required parameters are:

            mdot: accretion rate (g/cm)
                or
            lvisc: accretion luminosity in the disk due to viscous dissipation

            star: the central star

        '''

        # Disk Accretion
        if mdot is not None and lvisc is not None:
            raise Exception("Cannot specify both mdot and lvisc")
        elif mdot is not None:
            self.mdot = mdot
        elif lvisc is not None:
            self.lvisc = lvisc
        else:
            self.mdot = 0.

        # Central star
        self.star = star

        FlaredDisk.__init__(self, **kwargs)

    def __str__(self):
        string = FlaredDisk.__str__(self)
        string = string.replace('Flared', 'Alpha')
        string += " - Mdot: %.3e\n" % self.mdot
        string += " - Lvisc: %.3e\n" % self.lvisc
        return string

    def _check_all_set(self):

        FlaredDisk._check_all_set(self)

        if self.star is None:
            raise Exception("star is not set")

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

        int2 = integrate_powerlaw(self.rmin, self.rmax, 0.5 + self.p)
        int2 *= self.star.radius ** 0.5 * self.r_0 ** -self.p

        integral = (2. * pi) ** 1.5 * self.h_0 * (int1 - int2)

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

        # Geometrical factor
        rho *= (1. - np.sqrt(self.star.radius / grid.gw))

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

        print "Normalization factor for disk mass (should ideally be 1): %16.10f" % norm

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

        int2 = integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), -0.5 + self.p - self.beta)
        int2 *= self.star.radius ** 0.5 * self.r_0 ** (self.beta - self.p)

        return self.rho_0() * (int1 - int2)

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

        # Geometrical factor
        rho *= (1. - np.sqrt(self.star.radius / w))

        rho *= self.rho_0()

        # What about normalization

        return rho

    def vertical_cumulative_density(self, r, theta):

        density = self.vertical_profile(r, theta)

        d = r * np.radians(theta)

        tau = density * d

        tau[0] = 0.

        return tau

    def accretion_luminosity(self, grid):
        '''
        Return a luminosity grid for spherical polar coordinates

        Input is the position of the center of the grid cells in spherical
        polar coordinates (r, theta, phi), and the volume of the grid cells.
        '''

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(grid.shape)

        if 'lvisc' in self.__dict__ and self.lvisc == 0.:
            return np.zeros(grid.shape)

        if 'mdot' in self.__dict__ and self.mdot == 0.:
            return np.zeros(grid.shape)

        self._check_all_set()

        # Find disk scaleheight at each cylindrical radius
        h = self.h_0 * (grid.gw / self.r_0) ** self.beta

        # Find normalization constant
        if 'lvisc' in self.__dict__:

            int1 = integrate_powerlaw(self.rmin, self.rmax, -2.0)

            int2 = integrate_powerlaw(self.rmin, self.rmax, -2.5)
            int2 *= self.star.radius ** 0.5

            integral = (2. * pi) ** 1.5 * (int1 - int2)

            lvisc0 = self.lvisc / integral

        else:

            lvisc0 = 3. * G * self.star.mass * self.mdot \
                / np.sqrt(32. * pi ** 3.)

        # Find disk luminosity at all positions
        luminosity = lvisc0 / grid.gw ** 3 / h * grid.volumes \
                   * (1. - np.sqrt(self.star.radius / grid.gw)) \
                   * np.exp(-0.5 * (grid.gz / h) ** 2)

        # Truncate below rmin and above rmax
        if self.cylindrical_inner_rim:
            luminosity[grid.gw < self.rmin] = 0.
        else:
            luminosity[grid.gr < self.rmin] = 0.

        if self.cylindrical_outer_rim:
            luminosity[grid.gw > self.rmax] = 0.
        else:
            luminosity[grid.gr > self.rmax] = 0.

        print "Luminosity sum [actual] : %.3e" % np.sum(luminosity)
        print "Luminosity sum [theoretical] : %.3e" % self.lvisc

        return luminosity

    def __setattr__(self, attribute, value):

        if attribute == 'mdot':
            if 'lvisc' in self.__dict__:
                warnings.warn("Overriding value of lvisc with value derived from mdot")
                del self.lvisc
            object.__setattr__(self, attribute, value)
        elif attribute == 'lvisc':
            if 'mdot' in self.__dict__:
                warnings.warn("Overriding value of mdot with value derived from lvisc")
                del self.mdot
            object.__setattr__(self, attribute, value)
        elif attribute == 'dust' and value is not None:
            FreezableClass.__setattr__(self, 'dust', SphericalDust(value))
        else:
            FreezableClass.__setattr__(self, attribute, value)

    def __getattr__(self, attribute):

        if attribute == 'lvisc':
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute disk accretion luminosity")
            lvisc = G * self.star.mass * self.mdot / 2. \
                    * (3. / self.rmin - 3. / self.rmax \
                    - 2. * np.sqrt(self.star.radius / self.rmin ** 3.) \
                    + 2. * np.sqrt(self.star.radius / self.rmax ** 3.))
            return lvisc

        if attribute == 'mdot':
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute disk accretion rate")
            mdot = self.lvisc / G / self.star.mass * 2. \
                 / (3. / self.rmin - 3. / self.rmax \
                    - 2. * np.sqrt(self.star.radius / self.rmin ** 3.) \
                    + 2. * np.sqrt(self.star.radius / self.rmax ** 3.))
            return mdot

        raise AttributeError(attribute)
