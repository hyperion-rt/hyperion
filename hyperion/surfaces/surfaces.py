from __future__ import print_function

import os
import abc

import h5py
import numpy as np

from .surface_scattering_properties import SurfaceScatteringProperties

from ..util.functions import FreezableClass, is_numpy_array
from ..util.validator import validate_scalar


class Surface(FreezableClass):
    """
    This is the base class to represent surfaces.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, name=None, surface_properties=None):

        self.name = name
        self.surface_properties = surface_properties

    @property
    def name(self):
        """
        The name of the surface
        """
        return self._name

    @name.setter
    def name(self, value):
        if isinstance(value, basestring):
            self._name = value
        else:
            return TypeError("Surface name should be a string")

    @property
    def surface_properties(self):
        """
        Scattering/emission properties of the surface
        """
        return self._surface_properties

    @surface_properties.setter
    def surface_properties(self, value):
        if value is not None:
            if isinstance(value, SurfaceScatteringProperties) or isinstance(value, basestring):
                self._surface_properties = value
            else:
                raise TypeError("surface_properties should be a string or a SurfaceScatteringProperties instance")
        else:
            self._surface_properties = None

    def _check_all_set(self):
        if self.surface_properties is None:
            raise ValueError("surface_properties is not set")

    def write(self, handle):
        """
        Write the properties of the surface to an HDF5 group
        """

        self._check_all_set()

        # handle.attrs['surface_properties'] = self.surface_properties


class SphericalSurface(Surface):

    def __init__(self, name=None, surface_properties=None, radius=None,
                 position=(0., 0., 0.)):

        Surface.__init__(self, name=name,
                         surface_properties=surface_properties)

        self.radius = radius
        self.position = position

    @property
    def radius(self):
        """
        The radius of the source (cm)
        """
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            validate_scalar('radius', value, domain='positive')
        self._radius = value

    @property
    def position(self):
        '''
        The cartesian position of the surface ``(x, y, z)`` as a sequence of three floating-point values (cm)
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
        if self.radius is None:
            raise ValueError("radius is not set")
        Surface._check_all_set(self)

    def write(self, handle, name, copy=True, absolute_paths=False):
        """
        Write the properties of the surface to an HDF5 group

        TODO: allow dictionary of existing surfaces to be passed to avoid duplication. Also deal with linking better, like for dust.
        """

        self._check_all_set()

        # Write attributes
        g = handle.create_group(name)
        g.attrs['type'] = np.string_('sphere'.encode('utf-8'))
        g.attrs['x'] = self.position[0]
        g.attrs['y'] = self.position[1]
        g.attrs['z'] = self.position[2]
        g.attrs['r'] = self.radius
        Surface.write(self, g)

        # Write surface properties

        if copy:

            g_prop = g.create_group('surface_properties')

            if isinstance(self.surface_properties, basestring):
                self.surface_properties = SurfaceScatteringProperties(self.surface_properties)

            self.surface_properties.write(g_prop)

        else:

            if not isinstance(self.surface_properties, basestring):
                raise ValueError("cannot use copy=False unless surface_properties is set to a filename")

            if absolute_paths:
                path = os.path.abspath(self.surface_properties)
            else:
                # Relative path should be relative to input file, not current directory.
                path = os.path.relpath(self.surface_properties, os.path.dirname(g_prop.file))

            g['surface_properties'] = h5py.ExternalLink(path, '/')
