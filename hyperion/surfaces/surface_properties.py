from __future__ import print_statement

import atpy
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing
from ..version import __version__

# Bi-directional scattering distribution function

class SurfaceProperties(FreezableClass):

    def __init__(self):
        self.nu = None
        self.albedo = None
        self.theta_i = None
        self.theta_e = None
        self.g = None
        self.radiance = None

        self._freeze()
        
    @property
    def nu(self):
        """
        The frequencies for which the scattering properties are defined
        """
        return self._nu

    @nu.setter
    def nu(self, value):
        if value is None:
            self._nu = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("nu should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("nu should be monotonically increasing")
            if value[0] <= 0.:
                raise ValueError("nu should be strictly positive")
            self._nu = value

    @property
    def albedo(self):
        """
        The albedo of the surface
        """
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        if value is None:
            self._albedo = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("albedo should be a 1-D sequence")
            if np.any(value < 0.):
                raise ValueError("albedo should be positive")
            self._albedo = value

    @property
    def theta_i(self):
        """
        The incident theta angles for which the scattering properties are defined
        """
        return self._theta_i

    @theta_i.setter
    def theta_i(self, value):
        if value is None:
            self._theta_i = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("theta_i should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("theta_i should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("theta_i should be positive")
            self._theta_i = value

    @property
    def theta_e(self):
        """
        The emergent theta angles for which the scattering properties are defined
        """
        return self._theta_e

    @theta_e.setter
    def theta_e(self, value):
        if value is None:
            self._theta_e = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("theta_e should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("theta_e should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("theta_e should be positive")
            self._theta_e = value

    @property
    def g(self):
        """
        The emergent phase angles for which the scattering properties are defined
        """
        return self._g

    @g.setter
    def g(self, value):
        if value is None:
            self._g = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("g should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("g should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("g should be positive")
            self._g = value

    @property
    def radiance(self):
        """
        The bi-directional scattering distribution function.
        """
        return self._radiance

    @radiance.setter
    def radiance(self, value):
        if value is None:
            self._radiance = None
        else:
            if self.nu is None:
                raise Exception("nu has to be defined first")
            if self.theta_i is None:
                raise Exception("theta_i has to be defined first")
            if self.theta_e is None:
                raise Exception("theta_e has to be defined first")
            if self.g is None:
                raise Exception("g has to be defined first")
            if not is_numpy_array(value) or value.ndim != 4:
                raise ValueError("radiance should be a 4-d Numpy array")
            expected_shape = (len(self.nu), len(self.theta_i), len(self.theta_e), len(self.g))
            if value.shape != expected_shape:
                raise ValueError("radiance has an incorrect shape: {0:s} but expected {1:s}".format(value.shape, expected_shape))
            self._radiance = value

    def all_set(self):
        return self.nu is not None and \
               self.albedo is not None and \
               self.theta_i is not None and \
               self.theta_e is not None and \
               self.g is not None and \
               self.radiance is not None

    def write(self, handle, compression=True):

        if not self.all_set():
            raise Exception("Not all attributes of the surface properties are set")

        # Create dust table set
        ts = atpy.TableSet()

        # Add standard keywords to header
        ts.add_keyword('version', 1)
        ts.add_keyword('type', 1)
        ts.add_keyword('python_version', __version__)

        # Add optical properties
        topt = atpy.Table(name="optical_properties")
        topt.add_column('nu', self.nu)
        topt.add_column('albedo', self.albedo)
        ts.append(topt)

        # Add angles

        ti = atpy.Table(name="incident_angles")
        ti.add_column('i', self.theta_i)
        ts.append(ti)

        te = atpy.Table(name="emergent_angles")
        te.add_column('e', self.theta_e)
        ts.append(te)

        tg = atpy.Table(name="phase_angles")
        tg.add_column('g', self.g)
        ts.append(tg)

        # Write the table set to the HDF5 file
        ts.write(handle, type='hdf5')

        # Add the radiance as a dataset
        handle.create_dataset("radiance_pdf", data=self.radiance, compression=compression)
