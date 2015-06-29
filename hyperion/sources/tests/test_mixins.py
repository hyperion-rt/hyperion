from __future__ import print_function, division

import numpy as np
from astropy.tests.helper import pytest

from ..mixins import PositionMixin, VelocityMixin
from ..source import Source

class TestPositionMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, PositionMixin):
            pass

        self.p = SourceTest()

    def test_default(self):
        assert np.all(self.p.position == (0., 0., 0.))

    def test_set_non(self):
        self.p.position = None

    @pytest.mark.parametrize('position', [(0., 1., 2.),
                                          [3., 4., 5.],
                                          np.array([6., 7., 8.])])
    def test_set_valid(self, position):
        self.p.position = position

    @pytest.mark.parametrize('position', [(1., 2.),
                                          [3., 4., 5., 6],
                                          np.ones(5)])
    def test_set_invalid_len(self, position):
        with pytest.raises(ValueError) as exc:
            self.p.position = position
        assert exc.value.args[0] == "position should be a sequence of 3 values"

    def test_set_invalid_dim(self):
        with pytest.raises(ValueError) as exc:
            self.p.position = np.ones((3, 2))
        assert exc.value.args[0] == "position should be a 1-D sequence"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.p.position = ([0.], [1.,1.], [2.])
        assert exc.value.args[0] == "position should be a sequence of 3 scalar values"

    @pytest.mark.parametrize('position', ['a', 1., 4])
    def test_set_invalid_type(self, position):
        with pytest.raises(ValueError) as exc:
            self.p.position = position
        assert exc.value.args[0] == "position should be a tuple, list, or Numpy array"



class TestVelocityMixin(object):

    def setup_method(self, method):

        class SourceTest(Source, VelocityMixin):
            pass

        self.p = SourceTest()

    def test_default(self):
        assert self.p.velocity is None

    def test_set_non(self):
        self.p.velocity = None

    @pytest.mark.parametrize('velocity', [(0., 1., 2.),
                                          [3., 4., 5.],
                                          np.array([6., 7., 8.])])
    def test_set_valid(self, velocity):
        self.p.velocity = velocity

    @pytest.mark.parametrize('velocity', [(1., 2.),
                                          [3., 4., 5., 6],
                                          np.ones(5)])
    def test_set_invalid_len(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.p.velocity = velocity
        assert exc.value.args[0] == "velocity should be a sequence of 3 values"

    def test_set_invalid_dim(self):
        with pytest.raises(ValueError) as exc:
            self.p.velocity = np.ones((3, 2))
        assert exc.value.args[0] == "velocity should be a 1-D sequence"

    def test_set_invalid_shape(self):
        with pytest.raises(ValueError) as exc:
            self.p.velocity = ([0.], [1.,1.], [2.])
        assert exc.value.args[0] == "velocity should be a sequence of 3 scalar values"

    @pytest.mark.parametrize('velocity', ['a', 1., 4])
    def test_set_invalid_type(self, velocity):
        with pytest.raises(ValueError) as exc:
            self.p.velocity = velocity
        assert exc.value.args[0] == "velocity should be a tuple, list, or Numpy array"