from __future__ import print_function, division

import numpy as np
from astropy import log as logger

from ..util.functions import FreezableClass
from ..dust import SphericalDust
from ..util.validator import validate_scalar

from .core import Envelope, Density


class BipolarCavity(Density):
    r'''
    This class implements the density structure for a bipolar cavity
    in an envelope, with a density given by:

    .. math:: \rho(r) = \rho_0\,\left(\frac{r}{r_0}\right)^{-e} \\

    inside a volume defined by two parabolic surfaces with
    half-opening angle ``theta_0`` at ``r_0``.

    Once the :class:`~hyperion.densities.BipolarCavity` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> cavity = BipolarCavity()
        >>> cavity.theta_0 = 10.
        >>> cavity.power = 0.

    In most cases however, you should not have to instantiate a
    BipolarCavity class directly, but instead you should use the
    ``add_bipolar_cavity`` method on the Envelope classes (see for
    example :class:`UlrichEnvelope` or :class:`PowerLawEnvelope`
    classes).
    '''

    def __init__(self, theta_0=None, power=None, r_0=None,
                 rho_0=None, rho_exp=0.,
                 cap_to_envelope_density=False, dust=None):

        self.power = power
        self.theta_0 = theta_0
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.rho_exp = rho_exp

        self.cap_to_envelope_density = cap_to_envelope_density

        self._envelope = None

        self.dust = dust

        self._freeze()

    @property
    def theta_0(self):
        '''Cavity half-opening angle at ``r_0``'''
        return self._theta_0

    @theta_0.setter
    def theta_0(self, value):
        if value is not None:
            validate_scalar('theta_0', value, domain='positive')
        self._theta_0 = value

    @property
    def power(self):
        '''Power of the cavity shape'''
        return self._power

    @power.setter
    def power(self, value):
        if value is not None:
            validate_scalar('power', value, domain='real')
        self._power = value

    @property
    def r_0(self):
        '''radius at which ``theta_0`` and  ``rho_0`` are defined (cm)'''
        return self._r_0

    @r_0.setter
    def r_0(self, value):
        if value is not None:
            validate_scalar('r_0', value, domain='positive')
        self._r_0 = value

    @property
    def rho_0(self):
        '''density at ``r_0`` (g/cm^3)'''
        return self._rho_0

    @rho_0.setter
    def rho_0(self, value):
        if value is not None:
            validate_scalar('rho_0', value, domain='positive')
        self._rho_0 = value

    @property
    def rho_exp(self):
        '''density power-law exponent'''
        return self._rho_exp

    @rho_exp.setter
    def rho_exp(self, value):
        if value is not None:
            validate_scalar('rho_exp', value, domain='real')
        self._rho_exp = value

    @property
    def cap_to_envelope_density(self):
        '''
        Whether to always force the cavity density to be no larger
        than the density of the envelope.
        '''
        return self._cap_to_envelope_density

    @cap_to_envelope_density.setter
    def cap_to_envelope_density(self, value):
        if not isinstance(value, bool):
            raise ValueError("cap_to_envelope_density should be a boolean")
        self._cap_to_envelope_density = value

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

    def exists(self):
        return self.theta_0 > 0.

    def _check_all_set(self):

        if self.theta_0 is None:
            raise Exception("theta_0 is not set")
        if self.theta_0 > 0.:
            if self.power is None:
                raise Exception("power is not set")
            if self.r_0 is None:
                raise Exception("r_0 is not set")
            if self.rho_0 is None:
                raise Exception("rho_0 is not set")
            if self.rho_exp is None:
                raise Exception("rho_exp is not set")
            if self._envelope is None:
                raise Exception("envelope is not set")

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
            A 3-dimensional array containing the density of the
            bipolar cavity inside each cell. The shape of this array
            is the same as ``grid.shape``.
        '''

        self._check_all_set()

        if self.theta_0 == 0.:
            return np.zeros(grid.gr.shape)

        rho = self.rho_0 * np.abs(grid.gr / self.r_0) ** (-self.rho_exp)

        rho[grid.gr < self._envelope.rmin] = 0.
        rho[grid.gr > self._envelope.rmax] = 0.

        rho[self.mask(grid)] = 0.

        if self.cap_to_envelope_density:

            # Find envelope density
            envelope_density = self._envelope.density(grid, ignore_cavity=True)

            # Find cells where cavity density is larger than envelope and
            # reset
            reset = rho > envelope_density
            if np.all(reset):
                logger.warn("Bipolar cavity is denser than envelope everywhere, so will have no effect")
            rho[reset] = envelope_density[reset]

        return rho

    def mask(self, grid):
        '''
        Compute the shape of the bipolar cavity.

        Parameters
        ----------
        grid : :class:`~hyperion.grid.SphericalPolarGrid` or :class:`~hyperion.grid.CylindricalPolarGrid` instance.
            The spherical or cylindrical polar grid object containing
            information about the position of the grid cells.

        Returns
        -------
        mask : np.ndarray
            A 3-dimensional booleand array indicating whether we are
            inside or outside the bipolar cavity (True is inside). The
            shape of this array is the same as ``grid.shape``.
        '''

        if self.theta_0 == 0.:
            return np.ones(grid.gr.shape, dtype=bool)

        self._check_all_set()

        z0 = self.r_0 * np.cos(np.radians(self.theta_0))
        w0 = self.r_0 * np.sin(np.radians(self.theta_0))
        zcav = z0 * (grid.gw / w0) ** self.power

        mask = np.abs(grid.gz) < zcav

        return mask

    def __setattr__(self, attribute, value):

        if attribute == '_envelope':
            if value is not None and not isinstance(value, Envelope):
                raise ValueError("_envelope should be an instance of "
                                 "Envelope")

        FreezableClass.__setattr__(self, attribute, value)
