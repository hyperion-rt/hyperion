from __future__ import print_function, division

import numpy as np

from ..util.constants import pi
from ..util.functions import FreezableClass
from ..util.convenience import OptThinRadius
from ..util.integrate import integrate_powerlaw
from ..dust import SphericalDust
from ..util.logger import logger
from ..util.validator import validate_scalar


class FlaredDisk(FreezableClass):
    r'''
    This class implements the density structure for a flared axisymmatric
    disk, with a density given by:

    .. math:: \rho(R,z,\phi) = \rho_0^{\rm disk}\,\left(\frac{R_0}{R}\right)^{\beta - p}\,\exp{\left[-\frac{1}{2}\left(\frac{z}{h(R)}\right)^2\right]} \\

    Once the :class:`~hyperion.densities.FlaredDisk` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> from hyperion.util.constants import msun, au
        >>> disk = FlaredDisk()
        >>> disk.mass = 2. * msun
        >>> disk.rmin = 0.1 * au
        >>> disk.rmax = 100 * au
    '''

    def __init__(self, mass=0., rmin=None, rmax=None, p=-1,
                 beta=-1.25, h_0=None, r_0=None,
                 cylindrical_inner_rim=True,
                 cylindrical_outer_rim=True, dust=None):

        # Basic disk parameters
        self.mass = mass
        self.rmin = rmin
        self.rmax = rmax
        self.p = p
        self.beta = beta
        self.h_0 = h_0
        self.r_0 = r_0
        self.cylindrical_inner_rim = cylindrical_inner_rim
        self.cylindrical_outer_rim = cylindrical_outer_rim

        # Dust
        self.dust = dust

        self._freeze()

    @property
    def mass(self):
        '''total mass (g)'''
        return self._mass

    @mass.setter
    def mass(self, value):
        if value is not None:
            validate_scalar('mass', value, domain='positive')
        self._mass = value

    @property
    def rmin(self):
        '''inner radius (cm)'''
        return self._rmin

    @rmin.setter
    def rmin(self, value):
        if not isinstance(value, OptThinRadius) and value is not None:
            validate_scalar('rmin', value, domain='positive', extra=' or an OptThinRadius instance')
        self._rmin = value

    @property
    def rmax(self):
        '''outer radius (cm)'''
        return self._rmax

    @rmax.setter
    def rmax(self, value):
        if not isinstance(value, OptThinRadius) and value is not None:
            validate_scalar('rmax', value, domain='positive', extra=' or an OptThinRadius instance')
        self._rmax = value

    @property
    def p(self):
        '''surface density power-law exponent'''
        return self._p

    @p.setter
    def p(self, value):
        if value is not None:
            validate_scalar('p', value, domain='real')
        self._p = value

    @property
    def beta(self):
        '''scaleheight power-law exponent'''
        return self._beta

    @beta.setter
    def beta(self, value):
        if value is not None:
            validate_scalar('beta', value, domain='real')
        self._beta = value

    @property
    def h_0(self):
        '''scaleheight of the disk at ``r_0`` (cm)'''
        return self._h_0

    @h_0.setter
    def h_0(self, value):
        if value is not None:
            validate_scalar('h_0', value, domain='positive')
        self._h_0 = value

    @property
    def r_0(self):
        '''radius at which ``h_0`` is defined (cm)'''
        return self._r_0

    @r_0.setter
    def r_0(self, value):
        if value is not None:
            validate_scalar('r_0', value, domain='positive')
        self._r_0 = value

    @property
    def cylindrical_inner_rim(self):
        '''
        Whether the inner edge of the disk should be defined as a truncation
        in cylindrical or spherical polar coordinates
        '''
        return self._cylindrical_inner_rim

    @cylindrical_inner_rim.setter
    def cylindrical_inner_rim(self, value):
        if type(value) != bool:
            raise ValueError("cylindrical_inner_rim should be a boolean")
        self._cylindrical_inner_rim = value

    @property
    def cylindrical_outer_rim(self):
        '''
        Whether the outer edge of the disk should be defined as a truncation
        in cylindrical or spherical polar coordinates
        '''
        return self._cylindrical_outer_rim

    @cylindrical_outer_rim.setter
    def cylindrical_outer_rim(self, value):
        if type(value) != bool:
            raise ValueError("cylindrical_outer_rim should be a boolean")
        self._cylindrical_outer_rim = value

    @property
    def dust(self):
        '''dust properties (filename or dust object)'''
        return self._dust

    @dust.setter
    def dust(self, value):
        if isinstance(value, basestring):
            self._dust = SphericalDust(value)
        else:
            self._dust = value

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
            logger.warn("Ignoring disk, since rmax < rmin")
            return 0.

        int1 = integrate_powerlaw(self.rmin, self.rmax, 1.0 + self.p)
        int1 *= self.r_0 ** -self.p

        integral = (2. * pi) ** 1.5 * self.h_0 * int1

        return self.mass / integral

    def density(self, grid):
        '''
        Return the density grid

        Parameters
        ----------
        grid : :class:`~hyperion.grid.SphericalPolarGrid` or :class:`~hyperion.grid.CylindricalPolarGrid` instance.
            The spherical or cylindrical polar grid object containing
            information about the position of the grid cells.

        Returns
        -------
        rho : np.ndarray
            A 3-dimensional array containing the density of the disk inside
            each cell. The shape of this array is the same as
            ``grid.shape``.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            logger.warn("Ignoring disk, since rmax < rmin")
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

        if np.sum(rho * grid.volumes) == 0. and self.mass > 0:
            raise Exception("Discretized disk mass is zero, suggesting that the grid is too coarse")

        norm = self.mass / np.sum(rho * grid.volumes)

        logger.info("Disk density is being re-scaled by a factor of %.2f to give the correct mass." % norm)

        if norm > 1.1 or norm < 1. / 1.1:
            logger.warn("Re-scaling factor is significantly different from 1, which indicates that the grid may be too coarse to properly resolve the disk.")

        # Normalize to total disk mass
        rho = rho * norm

        return rho

    def midplane_cumulative_density(self, r):
        '''
        Find the cumulative column density as a function of radius.

        The cumulative density is measured outwards from the origin, and in
        the midplane.

        Parameters
        ----------
        r : np.ndarray
            Array of values of the radius up to which to tabulate the
            cumulative density.

        Returns
        -------
        rho : np.ndarray
            Array of values of the cumulative density.
        '''

        self._check_all_set()

        if self.rmax <= self.rmin:
            logger.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(r.shape)

        int1 = integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), self.p - self.beta)
        int1 *= self.r_0 ** (self.beta - self.p)

        return self.rho_0() * int1

    def _vertical_profile(self, r, theta):

        self._check_all_set()

        if self.rmax <= self.rmin:
            logger.warn("Ignoring disk, since rmax < rmin")
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
        '''
        Find the cumulative column density as a function of theta.

        Parameters
        ----------
        r : float
            The spherical radius at which to calculate the cumulative density.
        theta : np.ndarray
            The theta values at which to tabulate the cumulative density.

        Returns
        -------
        rho : np.ndarray
            Array of values of the cumulative density.
        '''
        density = self._vertical_profile(r, theta)

        d = r * np.radians(theta)

        tau = density * d

        tau[0] = 0.

        return tau
