import numpy as np

from hyperion.util.functions import B_nu, random_id
from hyperion.util.functions import FreezableClass
from hyperion.util.integrate import integrate_loglog

import atpy


class Source(FreezableClass):

    def __init__(self, luminosity=None, spectrum=None, temperature=None,
                 name=None):

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

        self._freeze()

    def check_all_set(self):
        if self.luminosity is None:
            raise Exception("luminosity is not set")

    def set_spectrum(self, spectrum):
        self.spectrum = spectrum

    def set_temperature(self, temperature):
        self.temperature = temperature

    def set_luminosity(self, luminosity):
        self.luminosity = luminosity

    def get_spectrum(self):
        self.check_all_set()
        if self.spectrum is not None:
            if isinstance(self.spectrum, atpy.Table):
                nu, fnu = self.spectrum.nu, self.spectrum.fnu
            elif type(self.spectrum) in [tuple, list]:
                nu, fnu = self.spectrum
            else:
                raise Exception("Spectrum should be tuple or ATpy table")
        elif self.temperature:
            nu = np.logspace(np.log10(3.e9), np.log10(3.e16))
            fnu = B_nu(nu, self.temperature)
        else:
            raise Exception("Not implemented")
        norm = integrate_loglog(nu, fnu)
        return nu, fnu / norm * self.luminosity

    def write_spectrum(self, handle):
        if self.spectrum:
            handle.attrs['spectrum'] = 'spectrum'
            if isinstance(self.spectrum, atpy.Table):
                self.spectrum.table_name = 'Spectrum'
                self.spectrum.write(handle, type='hdf5')
            else:
                table = atpy.Table(name='Spectrum')
                table.add_column('nu', self.spectrum[0])
                table.add_column('fnu', self.spectrum[1])
                table.write(handle, type='hdf5')
        elif self.temperature:
            handle.attrs['spectrum'] = 'temperature'
            handle.attrs['temperature'] = self.temperature
        else:
            handle.attrs['spectrum'] = 'lte'

    def has_lte_spectrum(self):
        return self.spectrum is None and self.temperature is None


class SpotSource(Source):

    def __init__(self, luminosity=None, longitude=None, latitude=None,
                 radius=None, spectrum=None, temperature=None, name=None):

        self.longitude = longitude
        self.latitude = latitude
        self.radius = radius
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.longitude is None:
            raise Exception("longitude is not set")
        if self.latitude is None:
            raise Exception("latitude is not set")
        if self.radius is None:
            raise Exception("radius is not set")
        if self.has_lte_spectrum():
            raise Exception("Spot source cannot have LTE spectrum")
        Source.check_all_set(self)

    def write(self, handle, name):
        self.check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = 'spot'
        g.attrs['luminosity'] = self.luminosity
        g.attrs['longitude'] = self.longitude
        g.attrs['latitude'] = self.latitude
        g.attrs['radius'] = self.radius
        self.write_spectrum(g)


class PointSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), spectrum=None,
                 temperature=None, name=None):
        self.position = position
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.has_lte_spectrum():
            raise Exception("Point source cannot have LTE spectrum")
        Source.check_all_set(self)

    def write(self, handle, name):
        self.check_all_set()
        g = handle.create_group(name)
        g.attrs['type'] = 'point'
        g.attrs['luminosity'] = self.luminosity
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        self.write_spectrum(g)


class SphericalSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), radius=None,
                 limb_darkening=False, spectrum=None, temperature=None,
                 name=None):
        self.position = position
        self.radius = radius
        self.limb = limb_darkening
        self.spots = []
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.radius is None:
            raise Exception("radius is not set")
        if self.limb is None:
            raise Exception("limb is not set")
        if self.has_lte_spectrum():
            raise Exception("Spherical source cannot have LTE spectrum")
        Source.check_all_set(self)

    def write(self, handle, name):

        self.check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'sphere'
        g.attrs['luminosity'] = self.luminosity
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        if self.limb:
            g.attrs['limb'] = 'yes'
        else:
            g.attrs['limb'] = 'no'
        self.write_spectrum(g)

        for i, spot in enumerate(self.spots):
            spot.write(g, 'Spot %i' % i)

    def add_spot(self, *args, **kwargs):
        self.spots.append(SpotSource(*args, **kwargs))


class ExternalSphericalSource(Source):

    def __init__(self, luminosity=None, position=(0., 0., 0.), radius=None,
                 spectrum=None, temperature=None, name=None):
        self.position = position
        self.radius = radius
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.position is None:
            raise Exception("position is not set")
        if self.radius is None:
            raise Exception("r is not set")
        if self.has_lte_spectrum():
            raise Exception("External spherical source cannot have LTE spectrum")
        Source.check_all_set(self)

    def write(self, handle, name):

        self.check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'extern_sph'
        g.attrs['luminosity'] = self.luminosity
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        self.write_spectrum(g)


class ExternalBoxSource(Source):

    def __init__(self, luminosity=None, bounds=None,
                 spectrum=None, temperature=None, name=None):
        self.bounds = bounds
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.bounds is None:
            raise Exception("bounds are not set")
        if self.has_lte_spectrum():
            raise Exception("External spherical source cannot have LTE spectrum")
        Source.check_all_set(self)

    def write(self, handle, name):

        self.check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'extern_box'
        g.attrs['luminosity'] = self.luminosity
        g.attrs['xmin'] = self.bounds[0][0]
        g.attrs['xmax'] = self.bounds[0][1]
        g.attrs['ymin'] = self.bounds[1][0]
        g.attrs['ymax'] = self.bounds[1][1]
        g.attrs['zmin'] = self.bounds[2][0]
        g.attrs['zmax'] = self.bounds[2][1]
        self.write_spectrum(g)


class MapSource(Source):

    def __init__(self, luminosity=None, map=None, spectrum=None,
                 temperature=None, name=None):
        self.map = map
        Source.__init__(self, luminosity=luminosity, spectrum=spectrum,
                        temperature=temperature, name=name)

    def check_all_set(self):
        if self.map is None:
            raise Exception("map is not set")
        if np.all(self.map == 0.):
            raise Exception("Luminosity map is zero everywhere")
        Source.check_all_set(self)

    def write(self, handle, name, grid, compression=True, map_dtype=float):

        self.check_all_set()

        g = handle.create_group(name)
        g.attrs['type'] = 'map'
        g.attrs['luminosity'] = self.luminosity
        grid.write_physical_array(g, self.map, "Luminosity map", dust=False,
                                  compression=compression,
                                  physics_dtype=map_dtype)
        self.write_spectrum(g)
