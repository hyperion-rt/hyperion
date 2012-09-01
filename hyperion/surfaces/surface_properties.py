from __future__ import print_function

import atpy
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing
from ..version import __version__


# The notation used here for the angles (mu0, mu, psi) is from Hapke (2012)

class SurfaceProperties(FreezableClass):

    def __init__(self):
        self.nu = None
        self.albedo = None
        self.mu0 = None
        self.mu = None
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
    def mu0(self):
        """
        The incident cos(theta) angles for which the scattering properties are defined
        """
        return self._mu0

    @mu0.setter
    def mu0(self, value):
        if value is None:
            self._mu0 = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("i should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("i should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("i should be positive")
            self._mu0 = value

    @property
    def mu(self):
        """
        The emergent cos(theta) angles for which the scattering properties are defined
        """
        return self._mu

    @mu.setter
    def mu(self, value):
        if value is None:
            self._mu = None
        else:
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("mu should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("mu should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("mu should be positive")
            self._mu = value

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
                raise ValueError("psi should be a 1-D sequence")
            if not monotonically_increasing(value):
                raise ValueError("psi should be monotonically increasing")
            if value[0] < 0.:
                raise ValueError("psi should be positive")
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
            if self.mu0 is None:
                raise Exception("mu0 has to be defined first")
            if self.mu is None:
                raise Exception("mu has to be defined first")
            if self.psi is None:
                raise Exception("psi has to be defined first")
            if not is_numpy_array(value) or value.ndim != 4:
                raise ValueError("radiance should be a 4-d Numpy array")
            expected_shape = (len(self.nu), len(self.mu0), len(self.mu), len(self.psi))
            if value.shape != expected_shape:
                raise ValueError("radiance has an incorrect shape: {0:s} but expected {1:s}".format(value.shape, expected_shape))
            self._radiance = value

    def all_set(self):
        return self.nu is not None and \
               self.albedo is not None and \
               self.mu0 is not None and \
               self.mu is not None and \
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

        tmu0 = atpy.Table(name="incident_angles")
        tmu0.add_column('mu0', self.mu0)
        ts.append(tmu0)

        tmu = atpy.Table(name="emergent_e_angles")
        tmu.add_column('mu', self.mu)
        ts.append(tmu)

        tpsi = atpy.Table(name="emergent_psi_angles")
        tpsi.add_column('psi', self.psi)
        ts.append(tpsi)

        # Write the table set to the HDF5 file
        ts.write(handle, type='hdf5')

        # Add the radiance as a dataset
        handle.create_dataset("radiance_pdf", data=self.radiance, compression=compression)
