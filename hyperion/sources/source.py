from __future__ import print_function, division

import atpy
import numpy as np

from ..grid.amr_grid import AMRGridView

from ..util.functions import B_nu, random_id, FreezableClass, \
                             is_numpy_array, bool2str, monotonically_increasing
from ..util.integrate import integrate_loglog
from ..util.validator import validate_scalar
from ..util.logger import logger


class Source(FreezableClass):
    '''
    This class is not meant to be used directly, but forms the basis for all other source types.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        if name:
            self.name = name
        else:
            self.name = random_id(length=8)

        self.peeloff = peeloff

        # Initialize attributes
        self.luminosity = None
        self.spectrum = None
        self.temperature = None

        # Prevent new attributes from being created
        self._freeze()

        # Set attributes from remaining keyword arguments
        for kwarg in kwargs:
            self.__setattr__(kwarg, kwargs[kwarg])

    @property
    def luminosity(self):
        '''
        The bolometric luminosity of the source (ergs/s)
        '''
        return self._luminosity

    @luminosity.setter
    def luminosity(self, value):
        if value is not None:
            validate_scalar('luminosity', value, domain='positive')
        self._luminosity = value

    @property
    def temperature(self):
        '''
        The temperature of the source (K)
        '''
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if value is not None:
            if hasattr(self, '_spectrum') and self._spectrum is not None:
                raise Exception("A spectrum has already been set, so cannot set a temperature")
            validate_scalar('temperature', value, domain='positive')
        self._temperature = value

    @property
    def spectrum(self):
        '''
        The spectrum of the source, specified either as an atpy.Table instance with ``'nu'`` and ``'fnu'`` columns, or as a ``(nu, fnu)`` tuple, where the frequency is given in Hz, and the flux is given as F_nu (units are unimportant since the normalization is set by the luminosity).
        '''
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):

        if value is not None:

            if hasattr(self, '_temperature') and self._temperature is not None:
                raise Exception("A temperature has already been set, so cannot set a spectrum")

            if isinstance(value, atpy.Table):

                if 'nu' not in value.columns:
                    raise TypeError("spectrum ATpy Table does not contain a"
                                    " 'nu' column")

                if 'fnu' not in value.columns:
                    raise TypeError("spectrum ATpy Table does not contain an"
                                    " 'fnu' column")

                nu, fnu = value['nu'], value['fnu']

            elif type(value) in (tuple, list):

                if len(value) == 2:
                    nu, fnu = value
                else:
                    raise TypeError("spectrum tuple or list should contain"
                                    " two elements")

                if type(nu) in [list, tuple]:
                    nu = np.array(nu, dtype=float)
                else:
                    nu = nu.astype(float)

                if type(fnu) in [list, tuple]:
                    fnu = np.array(fnu, dtype=float)
                else:
                    fnu = fnu.astype(float)

                if not is_numpy_array(nu) or nu.ndim != 1:
                    raise TypeError("nu should be a 1-D sequence")

                if not is_numpy_array(fnu) or fnu.ndim != 1:
                    raise TypeError("fnu should be a 1-D sequence")

                if nu.shape != fnu.shape:
                    raise TypeError("nu and fnu should have the same shape")

            else:

                raise TypeError('spectrum should be specified either as an '
                                'atpy.Table instance, or a tuple of two 1-D'
                                'Numpy arrays (nu, fnu) with the same length')

            # Check if frequency array has duplicate values
            if len(np.unique(nu)) != len(nu):
                raise ValueError("nu sequence contains duplicate values")

            # Check if spectrum needs sorting
            if not monotonically_increasing(nu):
                logger.warn("Spectrum is being re-sorted in order of increasing frequency")
                order = np.argsort(nu)
                nu = nu[order]
                fnu = fnu[order]

            self._spectrum = {'nu': nu, 'fnu': fnu}

        else:

            self._spectrum = value

    def _check_all_set(self):
        if self.luminosity is None:
            raise ValueError("luminosity is not set")

    def get_spectrum(self, nu_range=None):

        self._check_all_set()

        if self.spectrum is not None:
            nu, fnu = self.spectrum['nu'], self.spectrum['fnu']
            if nu_range is not None:
                raise NotImplemented("nu_range not yet implemented for spectrum")
        elif self.temperature is not None:
            if nu_range is None:
                raise ValueError("nu_range is needed for sources with Planck spectra")
            nu = np.logspace(np.log10(nu_range[0]), np.log10(nu_range[1]))
            nu[0] = nu_range[0]  # fix roundoff
            nu[-1] = nu_range[1]  # fix roundoff
            fnu = B_nu(nu, self.temperature)
        else:
            raise Exception("Not implemented")
        norm = integrate_loglog(nu, fnu)
        return nu, fnu / norm * self.luminosity

    def write(self, handle):

        self._check_all_set()

        handle.attrs['luminosity'] = self.luminosity

        handle.attrs['peeloff'] = np.string_(bool2str(self.peeloff))

        if self.spectrum is not None:
            handle.attrs['spectrum'] = np.string_('spectrum'.encode('utf-8'))
            table = atpy.Table(name='spectrum')
            table.add_column('nu', self.spectrum['nu'])
            table.add_column('fnu', self.spectrum['fnu'])
            table.write(handle, type='hdf5')
        elif self.temperature is not None:
            handle.attrs['spectrum'] = np.string_('temperature'.encode('utf-8'))
            handle.attrs['temperature'] = self.temperature
        else:
            handle.attrs['spectrum'] = np.string_('lte'.encode('utf-8'))

    def has_lte_spectrum(self):
        return self.spectrum is None and self.temperature is None


class SpotSource(Source):

    def __init__(self, name=None, peeloff=True, **kwargs):
        '''
        A spot on a spherical source

        Parameters
        ----------
        name : str, optional
            The name of the source
        peeloff : bool, optional
            Whether to peel-off photons from this source

        Any additional arguments are are used to initialize attributes.

        Attributes
        ----------
        luminosity : float
            The luminosity of the source (in ergs/s)
        spectrum : atpy.Table or tuple
            The spectrum of the source, specified either as:
                * an ATpy table with ``nu`` and ``fnu`` column
                * a ``(nu, fnu)`` tuple
            ``nu`` should be in Hz, and the units for ``fnu`` are not
            important, since the luminosity determined the absolute scaling``
        temperature : float
            The temperature of the source (in K)
        longitude : float
            The longitude of the spot on the spherical source (in degrees)
        latitude : float
            The latitude of the spot on the spherical source (in degrees)
        radius : float
            The radius of the spherical source (in cm)
        '''

        self.longitude = None
        self.latitude = None
        self.radius = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.longitude is None:
            raise ValueError("longitude is not set")
        if self.latitude is None:
            raise ValueError("latitude is not set")
        if self.radius is None:
            raise ValueError("radius is not set")
        if self.has_lte_spectrum():
            raise ValueError("Spot source cannot have LTE spectrum")

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = np.string_('spot'.encode('utf-8'))
        g.attrs['longitude'] = self.longitude
        g.attrs['latitude'] = self.latitude
        g.attrs['radius'] = self.radius
        Source.write(self, g)

    def __setattr__(self, attribute, value):

        if attribute == 'longitude' and value is not None:
            validate_scalar('longitude', value, domain=[0, 360])
        elif attribute == 'latitude' and value is not None:
            validate_scalar('latitude', value, domain=[-90, 90])
        elif attribute == 'radius' and value is not None:
            validate_scalar('radius', value, domain='positive')

        Source.__setattr__(self, attribute, value)


class PointSource(Source):
    '''
    A point source.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.position = (0., 0., 0.)

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def position(self):
        '''
        The cartesian position of the source ``(x, y, z)`` as a sequence of three floating-point values (cm)
        '''
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("position should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            else:
                raise ValueError("position should be a tuple, list, or Numpy array")
        self._position = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.position is None:
            raise ValueError("position is not set")
        if self.has_lte_spectrum():
            raise ValueError("Point source cannot have LTE spectrum")

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = np.string_('point'.encode('utf-8'))
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        Source.write(self, g)


class SphericalSource(Source):
    '''
    A spherical source

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.position = (0., 0., 0.)
        self.radius = None
        self.limb = False
        self._spots = []

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def radius(self):
        '''
        The radius of the source (cm)
        '''
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            validate_scalar('radius', value, domain='positive')
        self._radius = value

    @property
    def position(self):
        '''
        The cartesian position of the source ``(x, y, z)`` as a sequence of three floating-point values (cm)
        '''
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("position should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            else:
                raise ValueError("position should be a tuple, list, or Numpy array")
        self._position = value

    @property
    def limb(self):
        '''
        Whether to include limb darkening
        '''
        return self._limb

    @limb.setter
    def limb(self, value):
        if value is not None:
            if not type(value) == bool:
                raise ValueError("limb should be a boolean value (True/False)")
        self._limb = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.position is None:
            raise ValueError("position is not set")
        if self.radius is None:
            raise ValueError("radius is not set")
        if self.limb is None:
            raise ValueError("limb is not set")
        if self.has_lte_spectrum():
            raise ValueError("Spherical source cannot have LTE spectrum")

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = np.string_('sphere'.encode('utf-8'))
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        g.attrs['limb'] = np.string_(bool2str(self.limb))
        Source.write(self, g)

        for i, spot in enumerate(self._spots):
            spot.write(g, 'Spot %i' % i)

    def add_spot(self, *args, **kwargs):
        '''
        Add a spot to the source.

        All arguments are passed to :class:`~hyperion.sources.SpotSource`,
        so see that class for more details
        '''
        self._spots.append(SpotSource(*args, **kwargs))


class ExternalSphericalSource(Source):
    '''
    An spherical external source.

    This can be used for example to simulate the interstellar radiation
    field. This source is similar to :class:`~hyperion.sources.SphericalSource`
    but emits photons *inwards*.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.position = (0., 0., 0.)
        self.radius = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def radius(self):
        '''
        The radius of the source (cm)
        '''
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            validate_scalar('radius', value, domain='positive')
        self._radius = value

    @property
    def position(self):
        '''
        The cartesian position of the source ``(x, y, z)`` as a sequence of three floating-point values (cm)
        '''
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("position should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            else:
                raise ValueError("position should be a tuple, list, or Numpy array")
        self._position = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.position is None:
            raise ValueError("position is not set")
        if self.radius is None:
            raise ValueError("radius is not set")
        if self.has_lte_spectrum():
            raise ValueError("External spherical source cannot have LTE spectrum")

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = np.string_('extern_sph'.encode('utf-8'))
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        Source.write(self, g)


class ExternalBoxSource(Source):
    '''
    An rectangular external source.

    This can be used for example to simulate the interstellar radiation
    field. This source emits *inwards*.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.bounds = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def bounds(self):
        '''
        The cartesian boundaries of the rectangular box specified as
        ``[[xmin, xmax], [ymin, ymax], [zmin, zmax]]`` (cm)
        '''
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if np.array(value).shape != (3, 2):
                    raise ValueError("bounds should be a sequence of 3 pairs of values")
            elif is_numpy_array(value):
                if value.ndim != 2:
                    raise ValueError("bounds should be a 2-d array")
                if value.shape != (3, 2):
                    raise ValueError("bounds should have a shape of (3, 2)")
            else:
                raise ValueError("bounds should be a tuple, list, or Numpy array")
        self._bounds = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.bounds is None:
            raise ValueError("bounds is not set")
        if self.has_lte_spectrum():
            raise ValueError("External spherical source cannot have LTE spectrum")

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = np.string_('extern_box'.encode('utf-8'))
        g.attrs['xmin'] = self.bounds[0][0]
        g.attrs['xmax'] = self.bounds[0][1]
        g.attrs['ymin'] = self.bounds[1][0]
        g.attrs['ymax'] = self.bounds[1][1]
        g.attrs['zmin'] = self.bounds[2][0]
        g.attrs['zmax'] = self.bounds[2][1]
        Source.write(self, g)


class MapSource(Source):
    '''
    A diffuse source.

    This can be used for example to simulate the interstellar radiation
    field. This source emits *inwards*.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.map = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def map(self):
        '''
        The relative luminosity in each cell, given as a Numpy array or an AMRGridView instance
        '''
        return self._map

    @map.setter
    def map(self, value):
        if value is not None:
            if not is_numpy_array(value) and not isinstance(value, AMRGridView):
                raise ValueError("map should be a Numpy array or an AMRGridView instance")
        self._map = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.map is None:
            raise ValueError("map is not set")
        if is_numpy_array(self.map) and np.all(self.map == 0.):
            raise ValueError("map is zero everywhere")

    def write(self, handle, name, grid, compression=True, map_dtype=float):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = np.string_('map'.encode('utf-8'))
        grid.write_single_array(g, "Luminosity map", self.map,
                                  compression=compression,
                                  physics_dtype=map_dtype)
        Source.write(self, g)


class PlaneParallelSource(Source):
    '''
    A circular plane-parallel source.

    This source emits all photons in the same direction perpendicular to the
    plane of the source, and in one direction, like a beam.

    Parameters
    ----------
    name : str, optional
        The name of the source
    peeloff : bool, optional
        Whether to peel-off photons from this source

    Notes
    -----
    Any additional arguments are are used to initialize attributes.
    '''

    def __init__(self, name=None, peeloff=False, **kwargs):

        if peeloff:
            raise ValueError("Cannot peeloff plane parallel source")

        self.position = (0., 0., 0.)
        self.radius = None
        self.direction = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)

    @property
    def radius(self):
        '''
        The radius of the source (cm)
        '''
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            validate_scalar('radius', value, domain='positive')
        self._radius = value

    @property
    def position(self):
        '''
        The cartesian position of the source ``(x, y, z)`` as a sequence of three floating-point values (cm)
        '''
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("position should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            else:
                raise ValueError("position should be a tuple, list, or Numpy array")
        self._position = value

    @property
    def direction(self):
        '''
        The direction the photons should be emitted in ``(theta, phi)`` where ``theta`` and ``phi`` are spherical polar angles (rad)
        '''
        return self._direction

    @direction.setter
    def direction(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 2:
                    raise ValueError("direction should be a sequence of 2 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("direction should be a 1-D sequence")
                if len(value) != 2:
                    raise ValueError("direction should be a sequence of 2 values")
            else:
                raise ValueError("direction should be a tuple, list, or Numpy array")
        self._direction = value

    def _check_all_set(self):
        Source._check_all_set(self)
        if self.position is None:
            raise ValueError("position is not set")
        if self.radius is None:
            raise ValueError("radius is not set")
        if self.direction is None:
            raise ValueError("direction is not set")
        if self.has_lte_spectrum():
            raise ValueError("Point source cannot have LTE spectrum")

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = np.string_('plane_parallel'.encode('utf-8'))
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        g.attrs['theta'] = self.direction[0]
        g.attrs['phi'] = self.direction[1]
        Source.write(self, g)
