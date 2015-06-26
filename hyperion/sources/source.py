from __future__ import print_function, division

import numpy as np
from astropy.table import Table, Column

from ..grid.amr_grid import AMRGridView

from ..util.functions import B_nu, random_id, FreezableClass, \
    is_numpy_array, bool2str, str2bool, monotonically_increasing
from ..util.integrate import integrate_loglog
from ..util.validator import validate_scalar
from astropy import log as logger

from .mixins import (PositionMixin, VelocityMixin,
                     VectorPositionMixin, VectorVelocityMixin,
                     RadiusMixin)

def read_source(handle):
    source_type = handle.attrs['type'].decode('ascii')
    if source_type == 'spot':
        return SpotSource.read(handle)
    elif source_type == 'point':
        return PointSource.read(handle)
    elif source_type == 'point_collection':
        return PointSourceCollection.read(handle)
    elif source_type == 'sphere':
        return SphericalSource.read(handle)
    elif source_type == 'extern_sph':
        return ExternalSphericalSource.read(handle)
    elif source_type == 'extern_box':
        return ExternalBoxSource.read(handle)
    elif source_type == 'map':
        return MapSource.read(handle)
    elif source_type == 'plane_parallel':
        return PlaneParallelSource.read(handle)
    else:
        raise ValueError("Unexpected source type: {0}".format(source_type))


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

    _support_lte_spectrum = True

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

        # Hooks for mix-in classes
        self._read_hooks = []
        self._write_hooks = []
        self._required = ['luminosity']

        # Initialize mix-in classes
        import inspect
        for base_class in inspect.getmro(self.__class__):
            if 'Mixin' in base_class.__name__:
                base_class.__init__(self)

        # Prevent new attributes from being created
        self._freeze()

        # Set attributes from remaining keyword arguments
        for kwarg in kwargs:
            self.__setattr__(kwarg, kwargs[kwarg])

    @property
    def name(self):
        """
        A string identifying the source (optional)
        """
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

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

    def _read_luminosity(self, handle):
        self.luminosity = handle.attrs['luminosity']

    def _write_luminosity(self, handle):
        handle.attrs['luminosity'] = self.luminosity

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
        The spectrum of the source, specified either as an astropy.table.Table
        instance with ``'nu'`` and ``'fnu'`` columns, or as a ``(nu, fnu)``
        tuple, where the frequency is given in Hz, and the flux is given as
        F_nu (units are unimportant since the normalization is set by the
        luminosity).
        '''
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):

        if value is not None:

            if hasattr(self, '_temperature') and self._temperature is not None:
                raise Exception("A temperature has already been set, so cannot set a spectrum")

            if isinstance(value, Table):

                if 'nu' not in value.columns:
                    raise TypeError("spectrum Table does not contain a"
                                    " 'nu' column")

                if 'fnu' not in value.columns:
                    raise TypeError("spectrum Table does not contain an"
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
                                'astropy.table.Table instance, or a tuple '
                                'of two 1-D Numpy arrays (nu, fnu) with the '
                                'same length')

            # Check if frequency array has duplicate values
            if len(np.unique(nu)) != len(nu):
                raise ValueError("nu sequence contains duplicate values")

            # Check for any negative values

            if np.any(nu <= 0.):
                raise ValueError("nu should be strictly positive")

            if np.any(fnu < 0.):
                raise ValueError("fnu should be positive")

            # Check for any NaN or Inf values

            if np.any(np.isnan(nu) | np.isinf(nu)):
                raise ValueError("nu contains NaN/Inf values")

            if np.any(np.isnan(fnu) | np.isinf(fnu)):
                raise ValueError("fnu contains NaN/Inf values")

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

        if self.has_lte_spectrum() and not self._support_lte_spectrum:
            raise ValueError("{0} cannot have LTE spectrum".format(self.__class__.__name__))

        for attribute in self._required:
            if getattr(self, attribute) is None:
                raise ValueError("{0} is not set".format(attribute))

        if hasattr(self, '_check_all_set_specific'):
            self._check_all_set_specific()

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

    @classmethod
    def read(cls, handle):

        self = cls()

        self._read_luminosity(handle)

        if not handle.attrs['type'].decode('ascii') == self.short:
            raise ValueError("Source is not a {0}".format(self.__class__.__name__))

        self.name = handle.attrs['name'].decode('utf-8')

        self.peeloff = str2bool(handle.attrs['peeloff'])

        if handle.attrs['spectrum'] == b'spectrum':
            self.spectrum = Table(np.array(handle['spectrum']))
        elif handle.attrs['spectrum'] == b'temperature':
            self.temperature = handle.attrs['temperature']
        elif handle.attrs['spectrum'] == b'lte':
            pass
        else:
            raise ValueError('Unexpected value for `spectrum`: %s' % handle.attrs['spectrum'])

        # Read in attributes specific to the sub-class
        if hasattr(self, '_read_specific'):
            self._read_specific(handle)

        # Read in attributes from mix-ins
        for reader in self._read_hooks:
            reader(handle)

        return self

    def write(self, handle, name, grid=None, compression=True, dtype=float):

        self._check_all_set()

        g = handle.create_group(name)

        self._write_luminosity(g)

        g.attrs['type'] = np.string_('spot'.encode('utf-8'))

        g.attrs['name'] = np.string_(self.name.encode('utf-8'))

        g.attrs['peeloff'] = np.string_(bool2str(self.peeloff))

        if self.spectrum is not None:
            g.attrs['spectrum'] = np.string_('spectrum'.encode('utf-8'))
            table = Table()
            table.add_column(Column(data=self.spectrum['nu'], name='nu'))
            table.add_column(Column(data=self.spectrum['fnu'], name='fnu'))
            table.write(g, path='spectrum')
        elif self.temperature is not None:
            g.attrs['spectrum'] = np.string_('temperature'.encode('utf-8'))
            g.attrs['temperature'] = self.temperature
        else:
            g.attrs['spectrum'] = np.string_('lte'.encode('utf-8'))

        # Write out attributes specific to the sub-class
        if hasattr(self, '_write_specific'):
            if grid is None:
                self._write_specific(g)
            else:
                self._write_specific(g, grid=grid, compression=compression, dtype=dtype)

        # Write out attributes for mix-ins
        for writer in self._write_hooks:
            writer(g)

    def has_lte_spectrum(self):
        return self.spectrum is None and self.temperature is None


class SpotSource(Source, RadiusMixin):
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
    spectrum : astropy.table.Table or tuple
        The spectrum of the source, specified either as:
            * an Astropy Table with ``nu`` and ``fnu`` column
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

    label = 'spot'
    _support_lte_spectrum = False

    def __init__(self, name=None, peeloff=True, **kwargs):

        self.longitude = None
        self.latitude = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)
        self._required.extend(['longitude', 'latitude'])

    def _read_specific(self, handle):
        self.longitude = handle.attrs['longitude']
        self.latitude = handle.attrs['latitude']

    def _write_specific(self, handle):
        handle.attrs['longitude'] = self.longitude
        handle.attrs['latitude'] = self.latitude

    def __setattr__(self, attribute, value):

        if attribute == 'longitude' and value is not None:
            validate_scalar('longitude', value, domain=[0, 360])
        elif attribute == 'latitude' and value is not None:
            validate_scalar('latitude', value, domain=[-90, 90])

        Source.__setattr__(self, attribute, value)


class PointSource(Source, PositionMixin, VelocityMixin):
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

    label = 'point'
    _support_lte_spectrum = False


class PointSourceCollection(Source, VectorPositionMixin, VectorVelocityMixin):
    '''
    A point source collection.

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

    label = 'point_collection'
    _support_lte_spectrum = False

    @property
    def luminosity(self):
        '''
        The bolometric luminosity for the N sources as a 1-D Numpy array (ergs/s)
        '''
        return self._luminosity

    @luminosity.setter
    def luminosity(self, value):
        if value is not None:
            if is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("luminosity should be a 1-D array")
                if not np.all(value > 0.):
                    raise ValueError("luminosity should be positive")
                if self.position is not None and value.shape[0] != self.position.shape[0]:
                    raise ValueError("luminosity should be a 1-D array with the same number of rows as position")
            else:
                raise ValueError("luminosity should be a Numpy array")
        self._luminosity = value

    def _read_luminosity(self, handle):
        self.luminosity = np.array(handle['luminosity'])

    def _write_luminosity(self, handle):
        handle.create_dataset('luminosity', data=self.luminosity, compression=True)


class SphericalSource(Source, PositionMixin, RadiusMixin, VelocityMixin):
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

    short = 'sphere'
    _support_lte_spectrum = False

    def __init__(self, name=None, peeloff=True, **kwargs):
        self.limb = False
        self._spots = []
        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)
        self._required.append('limb')

    @property
    def limb(self):
        '''
        Whether to include limb darkening
        '''
        return self._limb

    @limb.setter
    def limb(self, value):
        if value is not None:
            if not isinstance(value, bool):
                raise ValueError("limb should be a boolean value (True/False)")
        self._limb = value

    def _read_specific(self, handle):

        self.limb = str2bool(handle.attrs['limb'])

        for group in handle:
            if 'Spot' in group:
                self._spots.append(SpotSource.read(handle[group]))

    def _write_specific(self, handle):

        handle.attrs['limb'] = np.string_(bool2str(self.limb))

        for i, spot in enumerate(self._spots):
            spot.write(handle, 'Spot %i' % i)

    def add_spot(self, *args, **kwargs):
        '''
        Add a spot to the source.

        All arguments are passed to :class:`~hyperion.sources.SpotSource`,
        so see that class for more details
        '''
        spot = SpotSource(*args, **kwargs)
        self._spots.append(spot)
        return spot


class ExternalSphericalSource(Source, PositionMixin, RadiusMixin):
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

    label = 'extern_sph'
    _support_lte_spectrum = False


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

    label = 'extern_box'
    _support_lte_spectrum = False

    def __init__(self, name=None, peeloff=True, **kwargs):
        self.bounds = None
        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)
        self._required.append('bounds')

    @property
    def bounds(self):
        '''
        The cartesian boundaries of the rectangular box specified
        as ``[[xmin, xmax], [ymin, ymax], [zmin, zmax]]`` (cm).
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

    def _read_specific(self, handle):
        self.bounds = [(handle.attrs['xmin'], handle.attrs['xmax']),
                       (handle.attrs['ymin'], handle.attrs['ymax']),
                       (handle.attrs['zmin'], handle.attrs['zmax'])]

    def _write_specific(self, handle):
        handle.attrs['xmin'] = self.bounds[0][0]
        handle.attrs['xmax'] = self.bounds[0][1]
        handle.attrs['ymin'] = self.bounds[1][0]
        handle.attrs['ymax'] = self.bounds[1][1]
        handle.attrs['zmin'] = self.bounds[2][0]
        handle.attrs['zmax'] = self.bounds[2][1]


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

    label = 'map'
    _support_lte_spectrum = True

    def __init__(self, name=None, peeloff=True, **kwargs):
        self.map = None
        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)
        self._required.append('map')

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

    def _check_all_set_specific(self):
        if is_numpy_array(self.map) and np.all(self.map == 0.):
            raise ValueError("map is zero everywhere")

    def _read_specific(self, handle):
        self.map = np.array(handle['Luminosity map'])

    def _write_specific(self, handle, grid=None, compression=True, dtype=float):
        grid.write_single_array(handle, "Luminosity map", self.map,
                                compression=compression,
                                physics_dtype=dtype)


class PlaneParallelSource(Source, PositionMixin, RadiusMixin):
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

    label = 'plane_parallel'
    _support_lte_spectrum = False

    def __init__(self, name=None, peeloff=False, **kwargs):

        if peeloff:
            raise ValueError("Cannot peeloff plane parallel source")

        self.direction = None

        Source.__init__(self, name=name, peeloff=peeloff, **kwargs)
        self._required.append('direction')

    @property
    def direction(self):
        '''
        The direction the photons should be emitted in ``(theta, phi)`` where
        ``theta`` and ``phi`` are spherical polar angles (rad)
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

    def _read_specific(self, handle):
        self.direction = (handle.attrs['theta'], handle.attrs['phi'])

    def _write_specific(self, handle):
        handle.attrs['theta'] = self.direction[0]
        handle.attrs['phi'] = self.direction[1]
