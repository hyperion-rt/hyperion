from __future__ import print_function, division

import numpy as np
from astropy.tests.helper import pytest

from ..mixins import PositionMixin, VelocityMixin, VectorPositionMixin, VectorVelocityMixin
from ..source import Source


class TestPositionMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, PositionMixin):
            pass

        self.s = SourceTest()

    def test_default(self):
        assert np.all(self.s.position == (0., 0., 0.))

    def test_set_non(self):
        self.s.position = None

    @pytest.mark.parametrize('position', [(0., 1., 2.),
                                          [3., 4., 5.],
                                          np.array([6., 7., 8.])])
    def test_set_valid(self, position):
        self.s.position = position

    @pytest.mark.parametrize('position', [(1., 2.),
                                          [3., 4., 5., 6],
                                          np.ones(5)])
    def test_set_invalid_len(self, position):
        with pytest.raises(ValueError) as exc:
            self.s.position = position
        assert exc.value.args[0] == "position should be a sequence of 3 values"

    def test_set_invalid_dim(self):
        with pytest.raises(ValueError) as exc:
            self.s.position = np.ones((3, 2))
        assert exc.value.args[0] == "position should be a 1-D sequence"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.s.position = ([0.], [1.,1.], [2.])
        assert exc.value.args[0] == "position should be a sequence of 3 scalar values"

    @pytest.mark.parametrize('position', ['a', 1., 4])
    def test_set_invalid_type(self, position):
        with pytest.raises(ValueError) as exc:
            self.s.position = position
        assert exc.value.args[0] == "position should be a tuple, list, or Numpy array"



class TestVelocityMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, VelocityMixin):
            pass

        self.s = SourceTest()

    def test_default(self):
        assert self.s.velocity is None

    def test_set_non(self):
        self.s.velocity = None

    @pytest.mark.parametrize('velocity', [(0., 1., 2.),
                                          [3., 4., 5.],
                                          np.array([6., 7., 8.])])
    def test_set_valid(self, velocity):
        self.s.velocity = velocity

    @pytest.mark.parametrize('velocity', [(1., 2.),
                                          [3., 4., 5., 6],
                                          np.ones(5)])
    def test_set_invalid_len(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = velocity
        assert exc.value.args[0] == "velocity should be a sequence of 3 values"

    def test_set_invalid_dim(self):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = np.ones((3, 2))
        assert exc.value.args[0] == "velocity should be a 1-D sequence"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = ([0.], [1.,1.], [2.])
        assert exc.value.args[0] == "velocity should be a sequence of 3 scalar values"

    @pytest.mark.parametrize('velocity', ['a', 1., 4])
    def test_set_invalid_type(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = velocity
        assert exc.value.args[0] == "velocity should be a tuple, list, or Numpy array"


class TestVectorPositionMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, VectorPositionMixin):
            pass

        self.s = SourceTest()

    def test_default(self):
        assert self.s.position is None

    def test_set_non(self):
        self.s.position = None

    def test_set_valid(self):
        self.s.position = np.ones((10,3))

    @pytest.mark.xfail
    def test_set_invalid_len(self):
        self.s.luminosity = np.ones(10)
        with pytest.raises(ValueError) as exc:
            self.s.position = np.ones((11, 3))
        assert exc.value.args[0] == "position should be a 2-D array with the same number of rows as luminosity"
        self.s.position = np.ones((10, 3))

    @pytest.mark.parametrize('position', [np.ones(3), np.ones((3,3,3))])
    def test_set_invalid_dim(self, position):
        with pytest.raises(ValueError) as exc:
            self.s.position = position
        assert exc.value.args[0] == "position should be a 2-D array"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.s.position = np.ones((3,4))
        assert exc.value.args[0] == "position should be an Nx3 array"

    @pytest.mark.parametrize('position', ['a', 1., 4, (1,2,3)])
    def test_set_invalid_type(self, position):
        with pytest.raises(ValueError) as exc:
            self.s.position = position
        assert exc.value.args[0] == "position should be a Numpy array"


class TestVectorVelocityMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, VectorVelocityMixin):
            pass

        self.s = SourceTest()

    def test_default(self):
        assert self.s.velocity is None

    def test_set_non(self):
        self.s.velocity = None

    def test_set_valid(self):
        self.s.velocity = np.ones((10,3))

    @pytest.mark.xfail
    def test_set_invalid_len(self):
        self.s.luminosity = np.ones(10)
        with pytest.raises(ValueError) as exc:
            self.s.velocity = np.ones((11, 3))
        assert exc.value.args[0] == "velocity should be a 2-D array with the same number of rows as luminosity"
        self.s.velocity = np.ones((10, 3))

    @pytest.mark.parametrize('velocity', [np.ones(3), np.ones((3,3,3))])
    def test_set_invalid_dim(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = velocity
        assert exc.value.args[0] == "velocity should be a 2-D array"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = np.ones((3,4))
        assert exc.value.args[0] == "velocity should be an Nx3 array"

    @pytest.mark.parametrize('velocity', ['a', 1., 4, (1,2,3)])
    def test_set_invalid_type(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.s.velocity = velocity
        assert exc.value.args[0] == "velocity should be a Numpy array"
