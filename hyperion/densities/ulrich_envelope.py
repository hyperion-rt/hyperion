from __future__ import print_function, division

import numpy as np
from astropy import log as logger

from ..dust import SphericalDust
from ..grid import SphericalPolarGrid
from ..util.constants import G, pi
from ..util.convenience import OptThinRadius
from ..util.integrate import integrate_powerlaw
from ..util.validator import validate_scalar

from .core import Envelope
from .bipolar_cavity import BipolarCavity


def delta_neg(r, q):

    rho = np.sqrt(-q.real ** 3)
    theta = np.arccos(r.real / rho)

    s = (rho ** (1. / 3.) * np.cos(theta / 3.)).astype(np.complex)
    s.imag = rho ** (1. / 3.) * np.sin(theta / 3.)

    t = (rho ** (1. / 3.) * np.cos(-theta / 3.)).astype(np.complex)
    t.imag = rho ** (1. / 3.) * np.sin(-theta / 3.)

    return s, t


def delta_pos(r, delta):

    dr = (r + np.sqrt(delta)).real
    s = dr ** (1. / 3.)
    neg = dr < 0.
    s[neg] = - (- dr[neg]) ** (1. / 3.)

    dr = (r - np.sqrt(delta)).real
    t = dr ** (1. / 3.)
    neg = dr < 0.
    t[neg] = - (- dr[neg]) ** (1. / 3.)

    return s, t


def cubic(c, d):
    '''
    Solve x**3 + c * x + d = 0
    '''

    c = c.astype(np.complex)
    d = d.astype(np.complex)

    q = c / 3.
    r = - d / 2.

    delta = q ** 3 + r ** 2

    pos = delta >= 0.

    s = np.zeros(c.shape, dtype=np.complex)
    t = np.zeros(c.shape, dtype=np.complex)

    if np.sum(pos) > 0:
        s[pos], t[pos] = delta_pos(r[pos], delta[pos])

    if np.sum(~pos) > 0:
        s[~pos], t[~pos] = delta_neg(r[~pos], q[~pos])

    x1 = s + t
    x2 = - (s + t) / 2. + np.sqrt(3.) / 2. * (s - t) * np.complex(0., 1.)
    x3 = - (s + t) / 2. - np.sqrt(3.) / 2. * (s - t) * np.complex(0., 1.)

    return x1, x2, x3


def solve_mu0(ratio, mu):

    x = cubic(ratio - 1., - mu * ratio)

    ambig = (x[0].imag == 0) & (x[1].imag == 0) & (x[2].imag == 0)

    v = x[0].real

    mask = ambig & (mu >= 0) & (x[0].real >= 0) \
                             & (x[1].real < 0) \
                             & (x[2].real < 0)
    v[mask] = x[0][mask].real
    mask = ambig & (mu >= 0) & (x[0].real < 0) \
                             & (x[1].real >= 0) \
                             & (x[2].real < 0)
    v[mask] = x[1][mask].real
    mask = ambig & (mu >= 0) & (x[0].real < 0) \
                             & (x[1].real < 0) \
                             & (x[2].real >= 0)
    v[mask] = x[2][mask].real

    mask = ambig & (mu < 0) & (x[0].real < 0) \
                            & (x[1].real >= 0) \
                            & (x[2].real >= 0)
    v[mask] = x[0][mask].real
    mask = ambig & (mu < 0) & (x[0].real >= 0) \
                            & (x[1].real < 0) \
                            & (x[2].real >= 0)
    v[mask] = x[1][mask].real
    mask = ambig & (mu < 0) & (x[0].real >= 0) \
                            & (x[1].real >= 0) \
                            & (x[2].real < 0)
    v[mask] = x[2][mask].real

    return v


class UlrichEnvelope(Envelope):
    r'''
    This class implements the density structure for a rotationally flattened
    and infalling envelope, with a density given by:

    .. math:: \rho(r,\theta) = \frac{\dot{M}_{\rm env}}{4\pi\left(G M_{\star} R_{\rm c}^3\right)^{1/2}}\left(\frac{r}{ R_{\rm c}}\right)^{-3/2}\left(1 + \frac{\mu}{\mu_0}\right)^{-1/2}\left(\frac{\mu}{\mu_0} + \frac{2\mu_0^2 R_{\rm c}}{r}\right)^{-1} = \rho_0^{\rm env}\left(\frac{r}{ R_{\rm c}}\right)^{-3/2}\left(1 + \frac{\mu}{\mu_0}\right)^{-1/2}\left(\frac{\mu}{\mu_0} + \frac{2\mu_0^2 R_{\rm c}}{r}\right)^{-1}

    where $\mu_0$ is given by the equation for the streamline:

    .. math:: \mu_0^3 + \mu_0\left(\frac{r}{ R_{\rm c}} - 1\right) - \mu\left(\frac{r}{ R_{\rm c}}\right) = 0

    Once the :class:`~hyperion.densities.UlrichEnvelope` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> from hyperion.util.constants import msun, au, pc
        >>> envelope = UlrichEnvelope()
        >>> envelope.rho_0 = 1.e-19
        >>> envelope.rmin = 0.1 * au
        >>> envelope.rmax = pc

    :class:`~hyperion.densities.UlrichEnvelope` instances can only be used
    with spherical polar grids at this time.
    '''

    def __init__(self, mdot=None, rho_0=None, rmin=None, rmax=None, rc=None,
                 ambient_density=0., star=None):

        # Basic envelope parameters
        self.rmin = rmin
        self.rmax = rmax
        self.rc = rc

        # Envelope Infall
        if mdot is not None and rho_0 is not None:
            raise Exception("Cannot specify both mdot and rho_0")
        else:
            self.mdot = mdot
            self.rho_0 = rho_0

        # Central star
        self.star = star

        # Cavity
        self.cavity = None

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
    def rc(self):
        '''inner radius (cm)'''
        return self._rc

    @rc.setter
    def rc(self, value):
        if value is not None:
            validate_scalar('rc', value, domain='positive')
        self._rc = value

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
    def mdot(self):
        '''infall rate (g/s)'''
        if self._mdot is not None:
            return self._mdot
        elif self._rho_0 is None:
            return None
        else:
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute infall rate")
            mdot = self.rho_0 * \
                   (4. * pi * np.sqrt(G * self.star.mass * self.rc ** 3.))
            return mdot

    @mdot.setter
    def mdot(self, value):
        if value is not None:
            validate_scalar('mdot', value, domain='positive')
            if self._rho_0 is not None:
                logger.warn("Overriding value of rho_0 with value derived from mdot")
                self._rho_0 = None
        self._mdot = value

    @property
    def rho_0(self):
        '''density factor (g/cm^3)'''
        if self._rho_0 is not None:
            return self._rho_0
        elif self._mdot is None:
            return None
        else:
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute density scaling")
            rho_0 = self.mdot / \
                    (4. * pi * np.sqrt(G * self.star.mass * self.rc ** 3.))
            return rho_0

    @rho_0.setter
    def rho_0(self, value):
        if value is not None:
            validate_scalar('rho_0', value, domain='positive')
            if self._mdot is not None:
                logger.warn("Overriding value of mdot with value derived from rho_0")
                self._mdot = None
        self._rho_0 = value

    @property
    def star(self):
        '''central star instance (needs a ``mass`` attribute)'''
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
            self._star = value

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

        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")
        if self.rc is None:
            raise Exception("rc is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner envelope radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Outer envelope radius needs to be computed first")

        if self.star is None:
            raise Exception("star is not set")

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
            logger.warn("Ignoring Ulrich envelope, since rmax < rmin")
            return np.zeros(grid.shape)

        # Find mu = cos(theta)
        mu = np.cos(grid.gt)

        # Find mu_0, the cosine of the angle of a streamline of infalling
        # particles at r=infinity.
        mu0 = solve_mu0(grid.gr / self.rc, mu)

        # Find Ulrich envelope density
        rho = self.rho_0 * (grid.gr / self.rc) ** -1.5 \
                * (1 + mu / mu0) ** -0.5 \
                * (mu / mu0 + 2. * mu0 ** 2 * self.rc / grid.gr) ** -1.

        mid1 = (np.abs(mu) < 1.e-10) & (grid.gr < self.rc)
        rho[mid1] = self.rho_0 / np.sqrt(grid.gr[mid1] / self.rc) \
                  / (1. - grid.gr[mid1] / self.rc) / 2.

        mid2 = (np.abs(mu) < 1.e-10) & (grid.gr > self.rc)
        rho[mid2] = self.rho_0 / np.sqrt(2. * grid.gr[mid2] / self.rc - 1) \
                  / (grid.gr[mid2] / self.rc - 1.)

        if np.any((np.abs(mu) < 1.e-10) & (grid.gr == self.rc)):
            raise Exception("Grid point too close to Ulrich singularity")

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        if not ignore_cavity and self.cavity is not None:
            mask = self.cavity.mask(grid)
            rho[~mask] = 0.

        return rho

    def outermost_radius(self, rho):
        '''
        Find the outermost radius at which the density of the envelope has
        fallen to `rho` (in the midplane).

        Parameters
        ----------
        rho : float
            The density for which to determine the radius

        Returns
        -------
        r : float
            The radius at which the density has fallen to ``rho``
        '''
        a, b, c, d = 2., -5., 4., -1 - (self.rho_0 / rho) ** 2
        p = (3 * a * c - b ** 2) / (3. * a ** 2)
        q = (2. * b ** 3 - 9. * a * b * c + 27. * a ** 2. * d) / (27 * a ** 3.)
        x1, _, _ = cubic(np.array([p]), np.array([q]))
        x1 = x1.real[0] - b / 3. / a
        return x1 * self.rc

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
            logger.warn("Ignoring Ulrich envelope, since rmax < rmin")
            return np.zeros(r.shape)

        gamma_0 = self.rmin / self.rc
        gamma_1 = r.clip(self.rmin, self.rmax) / self.rc

        rho = np.zeros(r.shape)

        if gamma_0 < 1.:

            rho[:] = self.rho_0 * self.rc \
                * (np.log((np.sqrt(gamma_1) + 1) / (1. - np.sqrt(gamma_1)))
                   - np.log((np.sqrt(gamma_0) + 1) / (1. - np.sqrt(gamma_0))))

            rho[gamma_1 >= 1.] = np.inf

        elif gamma_0 > 1:

            rho[:] = self.rho_0 * self.rc \
                * (np.log((np.sqrt(2. * gamma_1 - 1.) - 1.) / (np.sqrt(2. * gamma_1 - 1.) + 1.))
                   - np.log((np.sqrt(2. * gamma_0 - 1.) - 1.) / (np.sqrt(2. * gamma_0 - 1.) + 1.)))

        else:

            rho[:] = np.inf

        return rho

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
