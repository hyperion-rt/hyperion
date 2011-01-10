import numpy as np

from hyperion.util.functions import FreezableClass
from hyperion.dust import SphericalDust


class BipolarCavity(FreezableClass):

    def __init__(self, theta_0=0., power=1.5, r_0=None,
                 rho_0=None, rho_exp=0., rho_amb=0., envelope=None):
        '''
        Initialize a Bipolar Cavity instance. The required parameters are:

            theta_0: opening angle (degrees)
            power: shape exponent
            r_0: base radius in the cavity (cm)
            rho_0: base density in the cavity (g/cm^3)
            rho_exp: vertical density exponent
            rho_amb: ambient density (g/cm^3)

            envelope: the parent envelope

        '''

        self.power = power
        self.theta_0 = theta_0
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.rho_exp = rho_exp

        self.envelope = envelope

        self.dust = None

        self._freeze()

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
            if self.envelope is None:
                raise Exception("envelope is not set")

    def density(self, grid):
        '''
        Compute the density of the bipolar cavity
        '''

        self._check_all_set()

        if self.theta_0 == 0.:
            return np.zeros(grid.gr.shape)

        z0 = self.r_0 * np.cos(np.radians(self.theta_0))
        w0 = self.r_0 * np.sin(np.radians(self.theta_0))
        zcav = z0 * (grid.gw / w0) ** self.power

        rho = self.rho_0 * np.abs(grid.gr / self.r_0) ** (-self.rho_exp)
        if hasattr(self.envelope, 'rho_amb'):
            rho[rho < self.envelope.rho_amb] = self.envelope.rho_amb
        rho[np.abs(grid.gz) < zcav] = 0.

        rho[grid.gr < self.envelope.rmin] = 0.
        rho[grid.gr > self.envelope.rmax] = 0.

        rho[self.mask(grid)] = 0

        return rho

    def mask(self, grid):
        '''
        Compute the density of the bipolar cavity
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
        if attribute == 'dust' and value is not None and type(value) is str:
            FreezableClass.__setattr__(self, 'dust', SphericalDust(value))
        else:
            FreezableClass.__setattr__(self, attribute, value)
