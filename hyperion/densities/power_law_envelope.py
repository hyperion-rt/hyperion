from __future__ import print_function, division

import numpy as np
from astropy import log as logger

from ..dust import SphericalDust
from ..grid import SphericalPolarGrid
from ..util.constants import pi
from ..util.convenience import OptThinRadius
from ..util.integrate import integrate_powerlaw
from ..util.validator import validate_scalar

from .core import Envelope
from .bipolar_cavity import BipolarCavity


class PowerLawEnvelope(Envelope):
    r'''
    This class implements the density structure for a spherically symmetric
    power-law envelope, with a density given by:

    .. math:: \rho(r) = \rho_0^{\rm env}\,\left(\frac{r}{r_0}\right)^\gamma \\

    Once the :class:`~hyperion.densities.PowerLawEnvelope` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> from hyperion.util.constants import msun, au, pc
        >>> envelope = PowerLawEnvelope()
        >>> envelope.mass = msun
        >>> envelope.rmin = 0.1 * au
        >>> envelope.rmax = pc

    :class:`~hyperion.densities.PowerLawEnvelope` instances can only be used
    with spherical polar grids at this time.
    '''

    def __init__(self):

        # Basic envelope parameters
        self.rmin = None
        self.rmax = None
        self.power = None
        self.r_0 = None

        self.rho_0 = None
        self.mass = None

        # Cavity
        self.cavity = None

        # Central star
        self.star = None

        # Dust
        self.dust = None

        self._freeze()

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
    def power(self):
        '''density power-law exponent'''
        return self._power

    @power.setter
    def power(self, value):
        if value is not None:
            validate_scalar('power', value, domain='real')
        self._power = value

    @property
    def r_0(self):
        '''radius at which ``rho_0`` is defined (cm)'''
        return self._r_0

    @r_0.setter
    def r_0(self, value):
        if value is not None:
            validate_scalar('r_0', value, domain='positive')
        self._r_0 = value

    @property
    def mass(self):
        '''total mass (g)'''
        if self._mass is not None:
            return self._mass
        elif self._rho_0 is None:
            return None
        else:
            self._check_all_set()
            alpha = 3. + self.power
            mass = self.rho_0 / alpha * \
                (4. * pi * (self.rmax ** alpha - self.rmin ** alpha) / self.r_0 ** self.power)
            return mass

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
        '''density at ``r_0`` (g/cm^3)'''
        if self._rho_0 is not None:
            return self._rho_0
        elif self._mass is None:
            return None
        else:
            self._check_all_set()
            alpha = 3. + self.power
            rho_0 = self.mass * alpha / \
                (4. * pi * (self.rmax ** alpha - self.rmin ** alpha) / self.r_0 ** self.power)
            return rho_0

    @rho_0.setter
    def rho_0(self, value):
        if value is not None:
            validate_scalar('rho_0', value, domain='positive')
            if self._mass is not None:
                logger.warn("Overriding value of mass with value derived from rho_0")
                self._mass = None
        self._rho_0 = value

    @property
    def cavity(self):
        '''BipolarCavity instance'''
        return self._cavity

    @cavity.setter
    def cavity(self, value):
        if value is None:
            self._cavity = None
        else:
            if not isinstance(value, BipolarCavity):
                raise ValueError("cavity should be an instance of BipolarCavity")
            self._cavity = value
            self._cavity._envelope = self

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

    def density(self, grid, ignore_cavity=False):
        '''
        Return the density grid

        Parameters
        ----------
        grid : :class:`~hyperion.grid.SphericalPolarGrid` instance.
            The spherical polar grid object containing information about the
            position of the grid cells.

        Returns
        -------
        rho : np.ndarray
            A 3-dimensional array containing the density of the envelope
            inside each cell. The shape of this array is the same as
            ``grid.shape``.
        '''

        if not isinstance(grid, SphericalPolarGrid):
            raise TypeError("grid should be a SphericalPolarGrid instance")

        self._check_all_set()

        if self.rmax <= self.rmin:
            logger.warn("Ignoring power-law envelope, since rmax < rmin")
            return np.zeros(grid.shape)

        rho = self.rho_0 * (grid.gr / self.r_0) ** self.power

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        if self._rho_0 is None:
            norm = self.mass / np.sum(rho * grid.volumes)
            logger.info("Normalization factor for envelope mass: %5.2f" % norm)
            rho = rho * norm

        if not ignore_cavity and self.cavity is not None:
            mask = self.cavity.mask(grid)
            rho[~mask] = 0.

        return rho

    def outermost_radius(self, rho):
        '''
        Find the outermost radius at which the density of the envelope has
        fallen to `rho`.

        Parameters
        ----------
        rho : float
            The density for which to determine the radius

        Returns
        -------
        r : float
            The radius at which the density has fallen to ``rho``
        '''
        return self.r_0 * (rho / self.rho_0) ** (1. / self.power)

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
            logger.warn("Ignoring power-law envelope, since rmax < rmin")
            return np.zeros(r.shape)

        return self.rho_0 * integrate_powerlaw(self.rmin, r.clip(self.rmin, self.rmax), self.power) / self.r_0 ** self.power

    def add_bipolar_cavity(self):
        '''
        Add a bipolar cavity to the envelope.

        Returns
        -------
        cavity : :class:`BipolarCavity` instance
            The bipolar cavity instance, which can then be used to set the
            parameters of the cavity.
        '''

        if self.cavity is not None:
            raise Exception("Envelope already has a bipolar cavity")
        else:
            self.cavity = BipolarCavity()
            return self.cavity
