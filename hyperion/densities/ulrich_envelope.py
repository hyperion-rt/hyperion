import warnings

import numpy as np

from hyperion.util.constants import G, pi
from hyperion.util.functions import FreezableClass
from hyperion.densities.bipolar_cavity import BipolarCavity
from hyperion.util.convenience import OptThinRadius


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


class UlrichEnvelope(FreezableClass):

    def __init__(self, mdot=None, rho_0=None, rmin=None, rmax=None, rc=None,
                 ambient_density=0., star=None):
        '''
        Initialize a Ulrich envelope instance. The required parameters are:

            rmin: inner radius (cm)
            rmax: outer radius (cm)
            rc: centrifugal radius (cm)
            ambient_density: minimum density for the envelope (g/cm^3)
            mdot: infall rate (g/s)
            rho_0: density scaling (g/cm^3)
            star: the central star

        '''

        # Basic envelope parameters
        self.rmin = rmin
        self.rmax = rmax
        self.rc = rc
        self.rho_amb = ambient_density

        # Envelope Infall
        if mdot is not None and rho_0 is not None:
            raise Exception("Cannot specify both mdot and rho_0")

        if mdot is not None:
            self.mdot = mdot
        elif rho_0 is not None:
            self.rho_0 = rho_0
        else:
            self.mdot = 0.

        # Central star
        self.star = star

        # Cavity
        self.cavity = None

        # Dust
        self.dust = None

        self._freeze()

    def _check_all_set(self):

        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")
        if self.rc is None:
            raise Exception("rc is not set")
        if self.rho_amb is None:
            raise Exception("rho_amb is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner envelope radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Outer envelope radius needs to be computed first")

        if self.star is None:
            raise Exception("star is not set")

    def exists(self):
        return self.rho_0 > 0. or self.rho_amb > 0.

    def density(self, grid):
        '''
        Find the density of a Ulrich envelope

        Input is the position of the center of the grid cells in spherical
        polar coordinates (r, theta, phi), the volume of the grid cells, and a
        parameter dictionary
        '''

        self._check_all_set()

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

        rho[rho < self.rho_amb] = self.rho_amb

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        if self.cavity is not None:
            mask = self.cavity.mask(grid)
            rho[~mask] = 0.

        return rho

    def add_bipolar_cavity(self):
        if self.cavity is not None:
            raise Exception("Envelope already has a bipolar cavity")
        else:
            self.cavity = BipolarCavity()
            self.cavity.envelope = self
            return self.cavity

    def __setattr__(self, attribute, value):

        if attribute == 'mdot':
            if 'rho_0' in self.__dict__:
                warnings.warn("Overriding value of rho_0 with value derived from mdot")
                del self.rho_0
            object.__setattr__(self, attribute, value)
        elif attribute == 'rho_0':
            if 'mdot' in self.__dict__:
                warnings.warn("Overriding value of mdot with value derived from rho_0")
                del self.mdot
            object.__setattr__(self, attribute, value)
        else:
            FreezableClass.__setattr__(self, attribute, value)

    def __getattr__(self, attribute):

        if attribute == 'rho_0':
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute density scaling")
            rho_0 = self.mdot / \
                    (4. * pi * np.sqrt(G * self.star.mass * self.rc ** 3.))
            return rho_0

        if attribute == 'mdot':
            self._check_all_set()
            if self.star.mass is None:
                raise Exception("Stellar mass is undefined - cannot compute infall rate")
            mdot = self.rho_0 * \
                   (4. * pi * np.sqrt(G * self.star.mass * self.rc ** 3.))
            return mdot

        raise AttributeError(attribute)
