import numpy as np
from ..util.functions import is_numpy_array
from ..util.validator import validate_scalar


class PositionMixin(object):

    def __init__(self):
        self.position = (0., 0., 0.)
        self._read_hooks.append(self._read_position)
        self._write_hooks.append(self._write_position)
        self._required.append('position')

    @property
    def position(self):
        """
        The cartesian position of the source ``(x, y, z)`` as a sequence of three floating-point values (cm)
        """
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
                for v in value:
                    if not np.isscalar(v):
                        raise ValueError("position should be a sequence of 3 scalar values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("position should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("position should be a sequence of 3 values")
            else:
                raise ValueError("position should be a tuple, list, or Numpy array")
        self._position = value

    def _read_position(self, handle):
        self.position = (handle.attrs['x'], handle.attrs['y'], handle.attrs['z'])

    def _write_position(self, handle):
        handle.attrs['x'] = self.position[0]
        handle.attrs['y'] = self.position[1]
        handle.attrs['z'] = self.position[2]


class VelocityMixin(object):

    def __init__(self):
        self.velocity = None
        self._read_hooks.append(self._read_velocity)
        self._write_hooks.append(self._write_velocity)

    @property
    def velocity(self):
        """
        The cartesian velocity of the source ``(vx, vy, vz)`` as a sequence of three floating-point values (cm/s)
        """
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None:
            if type(value) in [tuple, list]:
                if len(value) != 3:
                    raise ValueError("velocity should be a sequence of 3 values")
            elif is_numpy_array(value):
                if value.ndim != 1:
                    raise ValueError("velocity should be a 1-D sequence")
                if len(value) != 3:
                    raise ValueError("velocity should be a sequence of 3 values")
            else:
                raise ValueError("velocity should be a tuple, list, or Numpy array")
        self._velocity = value

    def _read_velocity(self, handle):
        if 'vx' in handle.attrs:
            self.velocity = (handle.attrs['vx'], handle.attrs['vy'], handle.attrs['vz'])
        else:
            self.velocity = None

    def _write_velocity(self, handle):
        if self.velocity is not None:
            handle.attrs['vx'] = self.velocity[0]
            handle.attrs['vy'] = self.velocity[1]
            handle.attrs['vz'] = self.velocity[2]


class VectorPositionMixin(object):

    def __init__(self):
        self.position = None
        self._read_hooks.append(self._read_position)
        self._write_hooks.append(self._write_position)
        self._required.append('position')

    @property
    def position(self):
        """
        The cartesian position of the N sources ``(x, y, z)`` as a 2-D Numpy array with shape Nx3 (cm)
        """
        return self._position

    @position.setter
    def position(self, value):
        if value is not None:
            if is_numpy_array(value):
                if value.ndim != 2:
                    raise ValueError("position should be a 2-D array")
                if value.shape[1] != 3:
                    raise ValueError("position should be an Nx3 array")
                if self.luminosity is not None and value.shape[0] != self.luminosity.shape[0]:
                    raise ValueError("position should be a 2-D array with the same number of rows as luminosity")
            else:
                raise ValueError("position should be a Numpy array")
        self._position = value

    def _read_position(self, handle):
        self.position = np.array(handle['position'])

    def _write_position(self, handle):
        handle.create_dataset('position', data=self.position, compression=True)


class VectorVelocityMixin(object):

    def __init__(self):
        self.velocity = None
        self._read_hooks.append(self._read_velocity)
        self._write_hooks.append(self._write_velocity)

    @property
    def velocity(self):
        """
        The cartesian velocity of the source ``(vx, vy, vz)`` as a sequence of three floating-point values (cm/s)
        """
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None:
            if is_numpy_array(value):
                if value.ndim != 2:
                    raise ValueError("velocity should be a 2-D array")
                if value.shape[1] != 3:
                    raise ValueError("velocity should be a N x 3 array")
            else:
                raise ValueError("velocity should be a Numpy array")
        self._velocity = value

    def _read_velocity(self, handle):
        if 'velocity' in handle:
            self.velocity = handle['velocity']
        else:
            self.velocity = None

    def _write_velocity(self, handle):
        if self.velocity is not None:
            handle.create_dataset('velocity', data=self.velocity, compression=True)


class RadiusMixin(object):

    def __init__(self):
        self.radius = None
        self._read_hooks.append(self._read_radius)
        self._write_hooks.append(self._write_radius)
        self._required.append('radius')

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

    def _read_radius(self, handle):
        self.radius = handle.attrs['r']

    def _write_radius(self, handle):
        handle.attrs['r'] = self.radius
