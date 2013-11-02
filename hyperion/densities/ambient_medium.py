from __future__ import print_function, division

import numpy as np

from ..dust import SphericalDust
from ..grid import SphericalPolarGrid
from ..util.functions import FreezableClass
from ..util.convenience import OptThinRadius
from ..util.validator import validate_scalar

from .core import Density


class AmbientMedium(Density):
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

    By default, the ambient medium simply adds a constant density ``rho`` of
    dust to the whole model between the inner and outer radius. However, it
    is possible to pass components that should be subtracted from the
    constant density using the ``subtract=`` argument. In the following
    example::

        >>> e = PowerLawEnvelope()
        >>> AmbientMedium(subtract=[e])

    the ambient medium does not simply add a constant density ``rho`` of dust
    everywhere, but it adds dust such that the density never falls below
    ``rho`` between ``rmin`` and ``rmax`` - that is, it subtracts the density
    of component ``e`` from the ``rho``, with a minimum of zero. In regions
    where the density of component of ``e`` is larger than ``rho``, no dust is
    added.
    '''
    def __init__(self, rho=None, rmin=None, rmax=None, subtract=[]):

        # Basic ambient medium parameters
        self.rho = rho
        self.rmin = rmin
        self.rmax = rmax
        self.subtract = subtract

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
    def dust(self):
        '''dust properties (filename or dust object)'''
        return self._dust

    @property
    def subtract(self):
        '''Components to subtract from the ambient density medium'''
        return self._subtract

    @subtract.setter
    def subtract(self, value):
        if value is not None:
            if not isinstance(value, (list, tuple)):
                raise TypeError("subtract should be a list")
            for c in value:
                if not isinstance(c, Density):
                    raise TypeError("component in `subtract` should be a density instance")
                if isinstance(c, AmbientMedium):
                    raise TypeError("component in `subtract` cannot be an ambient density instance")

        self._subtract = value

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

        # Reset density outside rmin < r < rmax
        rho[grid.gr < self.rmin] = 0.
        rho[grid.gr > self.rmax] = 0.

        # Subtract specified components (if any)
        for component in self.subtract:
            rho -= component.density(grid)
        rho[rho < 0] = 0.

        return rho
