from __future__ import print_function

import atpy
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing
from ..version import __version__


# The notation used here for the angles (i, e, psi) is from Hapke (2012)

class SurfaceProperties(FreezableClass):

    def __init__(self):
        self.nu = None
        self.albedo = None
        self.i = None
        self.e = None
        self.psi = None
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
    def i(self):
        """
        The incident theta angles for which the scattering properties are defined
        """
        return self._i

    @i.setter
    def i(self, value):
        if value is None:
            self._i = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("i should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("i should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("i should be positive")
            self._i = value

    @property
    def e(self):
        """
        The emergent theta angles for which the scattering properties are defined
        """
        return self._e

    @e.setter
    def e(self, value):
        if value is None:
            self._e = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("e should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("e should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("e should be positive")
            self._e = value

    @property
    def psi(self):
        """
        The emergent psi for which the scattering properties are defined
        """
        return self._psi

    @psi.setter
    def psi(self, value):
        if value is None:
            self._psi = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("g should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("g should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("g should be positive")
            self._psi = value

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
            if self.i is None:
                raise Exception("i has to be defined first")
            if self.e is None:
                raise Exception("e has to be defined first")
            if self.psi is None:
                raise Exception("g has to be defined first")
            if not is_numpy_array(value) or value.ndim != 4:
                raise ValueError("radiance should be a 4-d Numpy array")
            expected_shape = (len(self.nu), len(self.i), len(self.e), len(self.psi))
            if value.shape != expected_shape:
                raise ValueError("radiance has an incorrect shape: {0:s} but expected {1:s}".format(value.shape, expected_shape))
            self._radiance = value

    def all_set(self):
        return self.nu is not None and \
               self.albedo is not None and \
               self.i is not None and \
               self.e is not None and \
               self.psi is not None and \
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
        ti.add_column('i', self.i)
        ts.append(ti)

        te = atpy.Table(name="emergent_e_angles")
        te.add_column('e', self.e)
        ts.append(te)

        tg = atpy.Table(name="emergent_psi_angles")
        tg.add_column('psi', self.psi)
        ts.append(tg)

        # Write the table set to the HDF5 file
        ts.write(handle, type='hdf5')

        # Add the radiance as a dataset
        handle.create_dataset("radiance_pdf", data=self.radiance, compression=compression)
