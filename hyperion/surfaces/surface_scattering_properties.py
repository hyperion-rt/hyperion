from __future__ import print_function

import h5py
import atpy
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, monotonically_increasing
from ..version import __version__


# The notation used here for the angles (mu0, mu, psi) is from Hapke (2012)

class SurfaceScatteringProperties(FreezableClass):

    def __init__(self, *args):
        self.nu = None
        self.albedo = None
        self.mu0 = None
        self.mu = None
        self.psi = None
        self.brdf = None

        self._freeze()

        if len(args) == 0:
            pass
        elif len(args) == 1:
            self.read(args[0])
        else:
            raise Exception("SphericalDust cannot take more than one argument")

    @property
    def nu(self):
        """
        Frequencies (in Hz).
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
        Surface albedo.
        """
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        if value is None:
            self._albedo = None
        else:
            if self.nu is None:
                raise Exception("nu has to be defined first")
            if type(value) in [list, tuple]:
                value = np.array(value)
            if not is_numpy_array(value) or value.ndim != 1:
                raise ValueError("albedo should be a 1-D sequence")
            if np.any(value < 0.):
                raise ValueError("albedo should be positive")
            if len(value) != len(self.nu):
                raise ValueError("albedo should have the same length as nu")
            self._albedo = value

    @property
    def mu0(self):
        """
        Incident mu_0 = cos(i) angles (in the [0:1] range).
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
        Emergent mu = cos(e) angles (in the [0:1] range).
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
        Emergent psi angles (in radians).
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
    def brdf(self):
        """
        The bidirectional reflectance distribution function.
        """
        return self._brdf

    @brdf.setter
    def brdf(self, value):
        if value is None:
            self._brdf = None
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
                raise ValueError("brdf should be a 4-d Numpy array")
            expected_shape = (len(self.nu), len(self.mu0), len(self.mu), len(self.psi))
            if value.shape != expected_shape:
                raise ValueError("brdf has an incorrect shape: {0:s} but expected {1:s}".format(value.shape, expected_shape))
            self._brdf = value

    def all_set(self):
        return self.nu is not None and \
               self.albedo is not None and \
               self.mu0 is not None and \
               self.mu is not None and \
               self.psi is not None and \
               self.brdf is not None

    def read(self, filename):
        """
        Read the properties from a file

        Parameters
        ----------
        filename : str
            The name of the file to read the properties from
        """

        ts = atpy.TableSet(filename)

        self.nu = ts['optical_properties']['nu']
        self.albedo = ts['optical_properties']['albedo']
        self.mu0 = ts['incident_angles']['mu0']
        self.mu = ts['emergent_e_angles']['mu']
        self.psi = ts['emergent_psi_angles']['psi']
        f = h5py.File(filename, 'r')
        self.brdf = np.array(f['brdf'])
        f.close()

    def write(self, filename_or_handle, compression=True):
        """
        Write the properties to a file or an HDF5 group

        Parameters
        ----------
        filename_or_handle : str or h5py.highlevel.File or h5py.highlevel.Group
            The filename of the file to create, or the HDF5 group to write the
            properties to.
        """

        if not self.all_set():
            raise Exception("Not all attributes of the surface properties are set")

        if isinstance(filename_or_handle, basestring):
            handle = h5py.File(filename_or_handle, 'w')
        else:
            handle = filename_or_handle

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

        # Add the brdf as a dataset
        handle.create_dataset("brdf", data=self.brdf, compression=compression)

        # If filename was specified, close the file
        if isinstance(filename_or_handle, basestring):
            handle.close()
