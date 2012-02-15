import numpy as np

from hyperion.util.functions import B_nu, random_id
from hyperion.util.functions import FreezableClass, is_numpy_array
from hyperion.util.integrate import integrate_loglog
from hyperion.util.constants import c

import atpy


class Source(FreezableClass):

    def __init__(self, luminosity=None, spectrum=None, temperature=None,
                 name=None, peeloff=True):

        self.luminosity = luminosity
        self.spectrum = None
        self.temperature = None

        if spectrum:
            self.set_spectrum(spectrum)
        elif temperature:
            self.set_temperature(temperature)

        if name:
            self.name = name
        else:
            self.name = random_id(length=8)

        self.peeloff = peeloff

        self._freeze()

    def _check_all_set(self):
        if self.luminosity is None:
            raise ValueError("luminosity is not set")

    def set_spectrum(self, spectrum):
        self.spectrum = spectrum

    def set_temperature(self, temperature):
        self.temperature = temperature

    def set_luminosity(self, luminosity):
        self.luminosity = luminosity

    def get_spectrum(self):

        self._check_all_set()

        if self.spectrum is not None:
            if isinstance(self.spectrum, atpy.Table):
                nu, fnu = self.spectrum.nu, self.spectrum.fnu
            elif type(self.spectrum) in [tuple, list]:
                nu, fnu = self.spectrum
            else:
                raise Exception("Spectrum should be tuple or ATpy table")
        elif self.temperature:
            nu = np.logspace(np.log10(c), np.log10(c * 1.e6))
            fnu = B_nu(nu, self.temperature)
        else:
            raise Exception("Not implemented")
        norm = integrate_loglog(nu, fnu)
        return nu, fnu / norm * self.luminosity

    def write(self, handle):

        self._check_all_set()

        handle.attrs['luminosity'] = self.luminosity

        handle.attrs['peeloff'] = 'yes' if self.peeloff else 'no'

        if self.spectrum is not None:
            handle.attrs['spectrum'] = 'spectrum'
            if isinstance(self.spectrum, atpy.Table):
                self.spectrum.table_name = 'Spectrum'
                self.spectrum.write(handle, type='hdf5')
            else:
                table = atpy.Table(name='Spectrum')
                table.add_column('nu', self.spectrum[0])
                table.add_column('fnu', self.spectrum[1])
                table.write(handle, type='hdf5')
        elif self.temperature is not None:
            handle.attrs['spectrum'] = 'temperature'
            handle.attrs['temperature'] = self.temperature
        else:
            handle.attrs['spectrum'] = 'lte'

    def has_lte_spectrum(self):
        return self.spectrum is None and self.temperature is None

    def __setattr__(self, attribute, value):

        if attribute == 'luminosity' and value is not None:

            if not np.isscalar(value):
                raise ValueError("luminosity should be a scalar value")
            if not np.isreal(value):
                raise ValueError("luminosity should be a numerical value")
            if value < 0.:
                raise ValueError("luminosity should be positive")

            object.__setattr__(self, attribute, value)

        elif attribute == 'spectrum' and value is not None:

            if isinstance(value, atpy.Table):

                if 'nu' not in value.columns:
                    raise TypeError("spectrum ATpy Table does not contain a"
                                    " 'nu' column")

                if 'fnu' not in value.columns:
                    raise TypeError("spectrum ATpy Table does not contain an"
                                    " 'fnu' column")

                object.__setattr__(self, attribute, value)

            elif type(value) in (tuple, list):

                if len(value) == 2:
                    nu, fnu = value
                else:
                    raise TypeError("spectrum tuple or list should contain"
                                    " two elements")

                if not is_numpy_array(nu) or nu.ndim != 1:
                    raise TypeError("nu should be specified as a 1-D Numpy"
                                    " array")

                if not is_numpy_array(fnu) or fnu.ndim != 1:
                    raise TypeError("fnu should be specified as a 1-D Numpy"
                                    " array")

                if nu.shape != fnu.shape:
                    raise TypeError("nu and fnu should have the same shape")

                # Reverse direction if needed
                if nu[-1] < nu[0]:
                    nu = nu[::-1]
                    fnu = fnu[::-1]

                object.__setattr__(self, attribute, (nu, fnu))

            else:

                raise TypeError('spectrum should be specified either as an '
                                'atpy.Table instance, or a tuple of two 1-D'
                                'Numpy arrays (nu, fnu) with the same length')

        else:

            object.__setattr__(self, attribute, value)

class SpotSource(Source):

    def __init__(self, luminosity=None, longitude=None, latitude=None,
                 radius=None, spectrum=None, temperature=None, name=None, peeloff=True):

        self.longitude = longitude
        self.latitude = latitude
        self.radius = radius
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.longitude is None:
            raise Exception("longitude is not set")
        if self.latitude is None:
            raise Exception("latitude is not set")
        if self.radius is None:
            raise Exception("radius is not set")
        if self.has_lte_spectrum():
            raise Exception("Spot source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = 'spot'
        g.attrs['longitude'] = self.longitude
        g.attrs['latitude'] = self.latitude
        g.attrs['radius'] = self.radius
        Source.write(self, g)


class PointSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), spectrum=None,
                 temperature=None, name=None, peeloff=True):
        self.position = position
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.position is None:
            raise ValueError("position is not set")
        if self.has_lte_spectrum():
            raise ValueError("Point source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = 'point'
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        Source.write(self, g)

    def __setattr__(self, attribute, value):

        if attribute == 'position' and value is not None:

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

        else:

            Source.__setattr__(self, attribute, value)

class SphericalSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), radius=None,
                 limb_darkening=False, spectrum=None, temperature=None,
                 name=None, peeloff=True):
        self.position = position
        self.radius = radius
        self.limb = limb_darkening
        self.spots = []
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.radius is None:
            raise Exception("radius is not set")
        if self.limb is None:
            raise Exception("limb is not set")
        if self.has_lte_spectrum():
            raise Exception("Spherical source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'sphere'
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        if self.limb:
            g.attrs['limb'] = 'yes'
        else:
            g.attrs['limb'] = 'no'
        Source.write(self, g)

        for i, spot in enumerate(self.spots):
            spot.write(g, 'Spot %i' % i)

    def add_spot(self, *args, **kwargs):
        self.spots.append(SpotSource(*args, **kwargs))


class ExternalSphericalSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), radius=None,
                 spectrum=None, temperature=None, name=None, peeloff=True):
        self.position = position
        self.radius = radius
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.radius is None:
            raise Exception("r is not set")
        if self.has_lte_spectrum():
            raise Exception("External spherical source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'extern_sph'
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        Source.write(self, g)


class ExternalBoxSource(Source):

    def __init__(self, luminosity=None, bounds=None,
                 spectrum=None, temperature=None, name=None, peeloff=True):
        self.bounds = bounds
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.bounds is None:
            raise Exception("bounds are not set")
        if self.has_lte_spectrum():
            raise Exception("External spherical source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'extern_box'
        g.attrs['xmin'] = self.bounds[0][0]
        g.attrs['xmax'] = self.bounds[0][1]
        g.attrs['ymin'] = self.bounds[1][0]
        g.attrs['ymax'] = self.bounds[1][1]
        g.attrs['zmin'] = self.bounds[2][0]
        g.attrs['zmax'] = self.bounds[2][1]
        Source.write(self, g)


class MapSource(Source):

    def __init__(self, luminosity=None, map=None, spectrum=None,
                 temperature=None, name=None, peeloff=True):
        self.map = map
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.map is None:
            raise Exception("map is not set")
        if np.all(self.map == 0.):
            raise Exception("Luminosity map is zero everywhere")
        Source._check_all_set(self)

    def write(self, handle, name, grid, compression=True, map_dtype=float):

        self._check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'map'
        grid.write_physical_array(g, self.map, "Luminosity map", dust=False,
                                  compression=compression,
                                  physics_dtype=map_dtype)
        Source.write(self, g)


class PlaneParallelSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), radius=None,
                 direction=(0., 0.), spectrum=None,
                 temperature=None, name=None, peeloff=False):
        if peeloff:
            raise Exception("Cannot peeloff plane parallel source")
        self.position = position
        self.radius = range
        self.direction = direction
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name, peeloff=peeloff)

    def _check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.radius is None:
            raise Exception("radius is not set")
        if self.direction is None:
            raise Exception("direction is not set")
        if self.has_lte_spectrum():
            raise Exception("Point source cannot have LTE spectrum")
        Source._check_all_set(self)

    def write(self, handle, name):
        self._check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = 'plane_parallel'
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        g.attrs['theta'] = self.direction[0]
        g.attrs['phi'] = self.direction[1]
        Source.write(self, g)
