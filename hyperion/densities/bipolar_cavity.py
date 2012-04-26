from __future__ import print_function, division

import numpy as np

from ..util.functions import FreezableClass
from ..dust import SphericalDust
from ..densities.envelope import Envelope
from ..util.validator import validate_scalar
from ..util.logger import logger

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
            cap_to_envelope_density: whether to use the envelope density as an upper limit

            envelope: the parent envelope

        '''

        self.power = power
        self.theta_0 = theta_0
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.rho_exp = rho_exp

        self.envelope = envelope

        self.dust = None

        self.cap_to_envelope_density = False

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

        rho = self.rho_0 * np.abs(grid.gr / self.r_0) ** (-self.rho_exp)

        rho[grid.gr < self.envelope.rmin] = 0.
        rho[grid.gr > self.envelope.rmax] = 0.

        rho[self.mask(grid)] = 0.

        if self.cap_to_envelope_density:

            # Find envelope density
            envelope_density = self.envelope.density(grid, ignore_cavity=True)

            # Find cells where cavity density is larger than envelope and
            # reset
            reset = rho > envelope_density
            if np.all(reset):
                logger.warn("Bipolar cavity is denser than envelope everywhere, so will have no effect")
            rho[reset] = envelope_density[reset]

        return rho

    def mask(self, grid):
        '''
        Compute the shape of the bipolar cavity
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

        if value is not None:

            # Dust specified as string
            if attribute == 'dust' and isinstance(value, basestring):
                FreezableClass.__setattr__(self, 'dust', SphericalDust(value))
                return

            # Positive scalars
            if attribute in ['theta_0', 'r_0', 'rho_0']:
                validate_scalar(attribute, value, domain='positive')

            # Scalars
            if attribute in ['power', 'rho_exp']:
                validate_scalar(attribute, value, domain='real')

            # Bipolar cavity
            if attribute == 'envelope':
                if not isinstance(value, Envelope):
                    raise ValueError("envelope should be an instance of "
                                     "Envelope")

        FreezableClass.__setattr__(self, attribute, value)
