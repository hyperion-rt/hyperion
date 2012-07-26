from __future__ import print_function, division

import numpy as np

from ..util.functions import FreezableClass
from ..util.convenience import OptThinRadius
from ..dust import SphericalDust
from ..util.validator import validate_scalar
from ..grid import SphericalPolarGrid


class AmbientMedium(FreezableClass):
    r'''
    This class implements the density structure for an ambient density
    medium defined by a constant density, and an inner and outer radius.

    Once the :class:`~hyperion.densities.AmbientMedium` class has been
    instantiated, the parameters for the density structure can be set via
    attributes::

        >>> from hyperion.util.constants import au, pc
        >>> ambient = AmbientMedium()
        >>> ambient.rho = 1.e-20  # cgs
        >>> ambient.rmin = 0.1 * au  # cm
        >>> ambient.rmax = pc  # cm

    :class:`~hyperion.densities.AmbientMedium` instances can only be used with
    spherical polar grids at this time.
    '''
    def __init__(self, rho=None, rmin=None, rmax=None):

        # Basic ambient medium parameters
        self.rho = rho
        self.rmin = rmin
        self.rmax = rmax

        # Dust
        self.dust = None

        self._freeze()

    @property
    def rho(self):
        '''Density of the ambient medium (g/cm^3)'''
        return self._rho

    @rho.setter
    def rho(self, value):
        if value is not None:
            validate_scalar('rho', value, domain='positive')
        self._rho = value

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

        if self.density is None:
            raise Exception("density is not set")
        if self.rmin is None:
            raise Exception("rmin is not set")
        if self.rmax is None:
            raise Exception("rmax is not set")

        if isinstance(self.rmin, OptThinRadius):
            raise Exception("Inner ambient medium radius needs to be computed first")
        if isinstance(self.rmax, OptThinRadius):
            raise Exception("Inner ambient medium radius needs to be computed first")

    def density(self, grid):
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

        rho = np.ones(grid.gr.shape) * self.rho

        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        return rho
