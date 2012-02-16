import warnings

import numpy as np

from hyperion.util.constants import pi
from hyperion.densities.envelope import Envelope
from hyperion.densities.bipolar_cavity import BipolarCavity
from hyperion.util.convenience import OptThinRadius
from hyperion.util.integrate import integrate_powerlaw
from hyperion.dust import SphericalDust


class PowerLawEnvelope(Envelope):

    def __init__(self, mass=None, rho_0=None, r_0=None, rmin=None, rmax=None, power=None, star=None):
        '''
        Initialize a power-law envelope instance. The required parameters are:

            mass: envelope mass (g)
            rmin: inner radius (cm)
            rmin: output radius (cm)
            power: power-law exponent

        '''

        # Basic envelope parameters
        self.rmin = rmin
        self.rmax = rmax
        self.power = power
        self.r_0 = r_0

        # Envelope Infall
        if mass is not None and rho_0 is not None:
            raise Exception("Cannot specify both mass and rho_0")

        if mass is not None:
            self.mass = mass
        elif rho_0 is not None:
            self.rho_0 = rho_0
        else:
            self.mass = 0.

        # Central star
        # self.star = star

        # Cavity
        self.cavity = None

        # Dust
        self.dust = None

        self._freeze()

    def _check_all_set(self):

        if self.r_0 is None:
            raise Exception("r_0 is not set")
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
        return self.rho_0 > 0.

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

        rho = self.rho_0 * (grid.gr/self.r_0) ** self.power

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        norm = self.mass / np.sum(rho * grid.volumes)

        print "Normalization factor for envelope mass: %5.2f" % norm

        rho = rho * norm

        if self.cavity is not None:
            mask = self.cavity.mask(grid)
            rho[~mask] = 0.

        return rho

    def outermost_radius(self, rho):
        '''
        Find the outermost radius at which the density of the envelope has
        fallen to `rho`.
        '''
        return self.r_0 * (rho / self.rho_0) ** (1. / self.power)

    def midplane_cumulative_density(self, r):
        '''
        Find the cumulative column density as a function of radius from the
        star in the midplane of a power-law envelope.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            warnings.warn("Ignoring power-law envelope, since rmax < rmin")
            return np.zeros(r.shape)

        return self.rho_0 * integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), self.power) / self.r_0 ** self.power

    def add_bipolar_cavity(self):
        if self.cavity is not None:
            raise Exception("Envelope already has a bipolar cavity")
        else:
            self.cavity = BipolarCavity()
            return self.cavity

    def __setattr__(self, attribute, value):

        if value is not None:

            # Dust specified as string
            if attribute == 'dust' and isinstance(value, basestring):
                Envelope.__setattr__(self, 'dust', SphericalDust(value))
                return

            # Positive scalars
            if attribute in ['r_0', 'rho_0', 'mass']:
                if not np.isscalar(value):
                    raise ValueError("{:s} should be a scalar value".format(attribute))
                if not np.isreal(value):
                    raise ValueError("{:s} should be a numerical value".format(attribute))
                if value < 0.:
                    raise ValueError("{:s} should be positive".format(attribute))

            # Scalars
            if attribute in ['power']:
                if not np.isscalar(value):
                    raise ValueError("{:s} should be a scalar value".format(attribute))
                if not np.isreal(value):
                    raise ValueError("{:s} should be a numerical value".format(attribute))

            # Radii (positive scalars or OptThinRadius instance)
            if attribute in ['rmin', 'rmax']:
                if not isinstance(value, OptThinRadius):
                    if not np.isscalar(value):
                        raise ValueError("{:s} should be a scalar value or an OptThinRadius instance".format(attribute))
                    if not np.isreal(value):
                        raise ValueError("{:s} should be a numerical value or an OptThinRadius instance".format(attribute))
                    if value < 0.:
                        raise ValueError("{:s} should be positive".format(attribute))

            # Bipolar cavity
            if attribute == 'cavity':
                if not isinstance(value, BipolarCavity):
                    raise ValueError("cavity should be an instance of BipolarCavity")

        if attribute == 'mass':
            if 'rho_0' in self.__dict__:
                warnings.warn("Overriding value of rho_0 with value derived from mass")
                del self.rho_0
            object.__setattr__(self, attribute, value)
        elif attribute == 'rho_0':
            if 'mass' in self.__dict__:
                warnings.warn("Overriding value of mass with value derived from rho_0")
                del self.mass
            object.__setattr__(self, attribute, value)
        else:
            Envelope.__setattr__(self, attribute, value)

        if attribute == 'cavity' and isinstance(value, BipolarCavity):
            self.cavity.envelope = self

    def __getattr__(self, attribute):

        if attribute == 'rho_0':
            self._check_all_set()
            alpha = 3. + self.power
            rho_0 = self.mass * alpha / \
                (4. * pi * (self.rmax ** alpha - self.rmin ** alpha) / self.r_0 ** self.power)
            return rho_0

        if attribute == 'mass':
            self._check_all_set()
            alpha = 3. + self.power
            mass = self.rho_0 / alpha * \
                (4. * pi * (self.rmax ** alpha - self.rmin ** alpha)  / self.r_0 ** self.power)
            return mass

        raise AttributeError(attribute)
