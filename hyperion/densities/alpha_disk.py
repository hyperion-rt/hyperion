from __future__ import print_function, division

import numpy as np
from astropy import log as logger

from ..dust import SphericalDust
from ..util.constants import pi, G
from ..util.convenience import OptThinRadius
from ..util.integrate import integrate_powerlaw
from ..util.validator import validate_scalar

from .core import Disk


class AlphaDisk(Disk):
    r'''
    This class implements the density structure for an alpha-accretion disk as
    implemented in `Whitney et al. (2003)
    <http://dx.doi.org/10.1086/375415>`_, with a density given by:

    .. math:: \rho(R,z,\phi) = \rho_0^{\rm disk}\,\left(1 - \sqrt{\frac{R_{\star}}{R}}\right)\left(\frac{R_0}{R}\right)^{\beta - p}\,\exp{\left[-\frac{1}{2}\left(\frac{z}{h(R)}\right)^2\right]} \\

    where

    .. math:: h(R) = h_0\left(\frac{R}{R_0}\right)^\beta

    The :math:`\rho_0^{\rm disk}` parameter does not need to be set directly
    (although it can be), and is instead automatically calculated when you set
    the disk mass. The exact equation relating :math:`\rho_0^{\rm disk}` to the
    disk mass can be found by integrating the equation for
    :math:`\rho(R,z,\phi)` over three dimensions and setting the result equal
    to the disk mass.

    Once the :class:`~hyperion.densities.AlphaDisk` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> from hyperion.util.constants import msun, au
        >>> disk = AlphaDisk()
        >>> disk.mass = 2. * msun
        >>> disk.rmin = 0.1 * au
        >>> disk.rmax = 100 * au

    The difference between :class:`~hyperion.densities.FlaredDisk` and
    :class:`~hyperion.densities.AlphaDisk` is that the latter includes an
    extra term in the density equation (:math:`1 - \sqrt{R_0/R}`)
    but most importantly that it allows for viscous accretion luminosity,
    specified either via an accretion rate, or an accretion luminosity. The
    relation between the accretion rate and the accretion luminosity in an
    infinitesimal volume is:

    .. math:: \frac{d\dot{E}_{\rm acc}}{dV} = \frac{3 G M_\star \dot{M}_{\rm acc}}{\sqrt{32 \pi^3} R^3 h(R)} \left(1 - \sqrt{\frac{R_{\star}}{R}}\right) \exp{\left[-\frac{1}{2}\left(\frac{z}{h(R)}\right)^2\right]}

    This is equation (4) from `Whitney et al. (2003)
    <http://dx.doi.org/10.1086/375415>`_. Once integrated over the whole disk,
    this gives a total luminosity of:

    .. math:: L_{\rm acc} = \frac{G\,M_\star\,M_{\rm acc}}{2} \left[3\left(\frac{1}{R_{\rm min}} - \frac{1}{R_{\rm max}}\right) - 2\left(\sqrt{\frac{R_\star}{R_{\rm min}^3}} - \sqrt{\frac{R_\star}{R_{\rm max}^3}}\right)\right]
    '''

    def __init__(self, mass=None, rho_0=None, rmin=None, rmax=None, p=-1,
                 beta=-1.25, h_0=None, r_0=None, cylindrical_inner_rim=True,
                 cylindrical_outer_rim=True, mdot=None, lvisc=None, star=None,
                 dust=None):

        # Start off by initializing mass and rho_0
        self.mass = None
        self.rho_0 = None

        # Basic disk parameters
        self.rmin = rmin
        self.rmax = rmax
        self.p = p
        self.beta = beta
        self.h_0 = h_0
        self.r_0 = r_0
        self.cylindrical_inner_rim = cylindrical_inner_rim
        self.cylindrical_outer_rim = cylindrical_outer_rim

        # Disk mass
        if mass is not None and rho_0 is not None:
            raise Exception("Cannot specify both mass and rho_0")
        elif mass is not None:
            self.mass = mass
        elif rho_0 is not None:
            self.rho_0 = rho_0

        # Disk Accretion
        if mdot is not None and lvisc is not None:
            raise Exception("Cannot specify both mdot and lvisc")
        self.mdot = mdot
        self.lvisc = lvisc

        # Central star
        self.star = star

        # Dust
        self.dust = dust

        self._freeze()

    @property
    def mass(self):
        """
        Total disk mass (g)
        """
        if self._mass is not None:
            return self._mass
        elif self._rho_0 is None:
            return None
        else:

            self._check_all_set()

            if self.rmax <= self.rmin:
                return 0.

            int1 = integrate_powerlaw(self.rmin, self.rmax, 1.0 + self.p)
            int1 *= self.r_0 ** -self.p

            int2 = integrate_powerlaw(self.rmin, self.rmax, 0.5 + self.p)
            int2 *= self.star.radius ** 0.5 * self.r_0 ** -self.p

            integral = (2. * pi) ** 1.5 * self.h_0 * (int1 - int2)

            return self._rho_0 * integral

    @mass.setter
    def mass(self, value):
        if value is not None:
            validate_scalar('mass', value, domain='positive')
            if self._rho_0 is not None:
                logger.warn("Overriding value of rho_0 with value derived from mass")
                self._rho_0 = None
        self._mass = value

    @property
    def rho_0(self):
        """
        Scale-factor for the disk density (g/cm^3)
        """
        if self._rho_0 is not None:
            return self._rho_0
        elif self._mass is None:
            return None
        else:

            self._check_all_set()

            if self.rmax <= self.rmin:
                return 0.

            int1 = integrate_powerlaw(self.rmin, self.rmax, 1.0 + self.p)
            int1 *= self.r_0 ** -self.p

            int2 = integrate_powerlaw(self.rmin, self.rmax, 0.5 + self.p)
            int2 *= self.star.radius ** 0.5 * self.r_0 ** -self.p

            integral = (2. * pi) ** 1.5 * self.h_0 * (int1 - int2)

            return self._mass / integral

    @rho_0.setter
    def rho_0(self, value):
        if value is not None:
            validate_scalar('rho_0', value, domain='positive')
            if self._mass is not None:
                logger.warn("Overriding value of mass with value derived from rho_0")
                self._mass = None
        self._rho_0 = value

    @property
    def rmin(self):
        '''inner radius (cm)'''
        if isinstance(self._rmin, OptThinRadius):
            return self._rmin.evaluate(self.star, self.dust)
        else:
            return self._rmin

    @rmin.setter
    def rmin(self, value):
        if not isinstance(value, OptThinRadius) and value is not None:
            validate_scalar('rmin', value, domain='positive', extra=' or an OptThinRadius instance')
        self._rmin = value

    @property
    def rmax(self):
        '''outer radius (cm)'''
        if isinstance(self._rmax, OptThinRadius):
            return self._rmax.evaluate(self.star, self.dust)
        else:
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
        if not isinstance(value, bool):
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
        if not isinstance(value, bool):
            raise ValueError("cylindrical_outer_rim should be a boolean")
        self._cylindrical_outer_rim = value

    @property
    def mdot(self):
        '''accretion rate (g/s)'''
        if self._mdot is not None:
            return self._mdot
        elif self._lvisc is None:
            return None
        else:
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute disk accretion rate")
            mdot = self.lvisc / G / self.star.mass * 2. \
                 / (3. / self.rmin - 3. / self.rmax \
                    - 2. * np.sqrt(self.star.radius / self.rmin ** 3.) \
                    + 2. * np.sqrt(self.star.radius / self.rmax ** 3.))
            return mdot

    @mdot.setter
    def mdot(self, value):
        if value is not None:
            validate_scalar('mdot', value, domain='positive')
            if self._lvisc is not None:
                logger.warn("Overriding value of lvisc with value derived from mdot")
                self._lvisc = None
        self._mdot = value

    @property
    def lvisc(self):
        '''viscous accretion luminosity (ergs/s)'''
        if self._lvisc is not None:
            return self._lvisc
        elif self._mdot is None:
            return None
        else:
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute disk accretion luminosity")
            lvisc = G * self.star.mass * self.mdot / 2. \
                    * (3. / self.rmin - 3. / self.rmax \
                    - 2. * np.sqrt(self.star.radius / self.rmin ** 3.) \
                    + 2. * np.sqrt(self.star.radius / self.rmax ** 3.))
            return lvisc

    @lvisc.setter
    def lvisc(self, value):
        if value is not None:
            validate_scalar('lvisc', value, domain='positive')
            if self._mdot is not None:
                logger.warn("Overriding value of mdot with value derived from lvisc")
                self._mdot = None
        self._lvisc = value

    @property
    def star(self):
        '''central star instance (needs ``mass`` and ``radius`` attributes)'''
        return self._star

    @star.setter
    def star(self, value):
        if value is None:
            self._star = None
        else:
            try:
                value.mass
            except AttributeError:
                raise ValueError("star should have a ``mass`` attribute")
            try:
                value.radius
            except AttributeError:
                raise ValueError("star should have a ``radius`` attribute")
            self._star = value

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
        string = "= Alpha disk =\n"
        string += " - M_disk: %.3e\n" % self.mass
        string += " - R_min: %.3e\n" % self.rmin
        string += " - R_min: %.3e\n" % self.rmax
        string += " - p: %.3f\n" % self.p
        string += " - beta: %.3f\n" % self.beta
        string += " - h_0: %.3e\n" % self.h_0
        string += " - r_0: %.3e\n" % self.r_0
        string += " - Mdot: %.3e\n" % self.mdot
        string += " - Lvisc: %.3e\n" % self.lvisc
        return string

    def _check_all_set(self):

        if self._mass is None and self._rho_0 is None:
            raise Exception("either mass or rho_0 should be set")

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

        if self.star is None:
            raise Exception("star is not set")

        if self._lvisc is None and self._mdot is None:
            raise Exception("either lvisc or mdot should be set")

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
        rho *= self.rho_0

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

        int2 = integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), -0.5 + self.p - self.beta)
        int2 *= self.star.radius ** 0.5 * self.r_0 ** (self.beta - self.p)

        return self.rho_0 * (int1 - int2)

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

        # Geometrical factor
        rho *= (1. - np.sqrt(self.star.radius / w))

        rho *= self.rho_0

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

    def accretion_luminosity(self, grid):
        '''
        Return the viscous accretion luminosity grid

        Parameters
        ----------
        grid : :class:`~hyperion.grid.SphericalPolarGrid` or :class:`~hyperion.grid.CylindricalPolarGrid` instance.
            The spherical or cylindrical polar grid object containing
            information about the position of the grid cells.

        Returns
        -------
        lvisc : np.ndarray
            A 3-dimensional array containing the viscous accretion luminosity
            of the disk inside each cell. The shape of this array is the same
            as ``grid.shape``.
        '''

        if self.rmax <= self.rmin:
            logger.warn("Ignoring disk, since rmax < rmin")
            return np.zeros(grid.shape)

        if self.lvisc == 0.:
            return np.zeros(grid.shape)

        if self.mdot == 0.:
            return np.zeros(grid.shape)

        self._check_all_set()

        # Find disk scaleheight at each cylindrical radius
        h = self.h_0 * (grid.gw / self.r_0) ** self.beta

        # Find normalization constant
        if self.lvisc is not None:

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

        logger.info("Luminosity sum [actual] : %.3e" % np.sum(luminosity))
        logger.info("Luminosity sum [theoretical] : %.3e" % self.lvisc)

        return luminosity

    def scale_height_at(self, r):
        '''
        Return the scaleheight of the disk at radius `r`
        '''
        return self.h_0 * (r / self.r_0) ** self.beta
