from __future__ import print_function, division

import pytest
import numpy as np

from ..optical_properties import OpticalProperties


def test_init():
    OpticalProperties()

VECTOR_ATTRIBUTES = ['nu', 'chi', 'albedo', 'mu']
ARRAY_ATTRIBUTES = ['P1', 'P2', 'P3', 'P4']


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_list(attribute):
    o = OpticalProperties()
    setattr(o, attribute, [0.1, 0.2, 0.3])


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_array(attribute):
    o = OpticalProperties()
    setattr(o, attribute, np.array([0.1, 0.2, 0.3]))


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_type1(attribute):
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, 'hello')
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_type2(attribute):
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, 0.5)
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_shape1(attribute):
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, [[0., 1.], [0.5, 1.]])
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_shape2(attribute):
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, np.array([[0., 1.], [0.5, 1.]]))
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_order(attribute):
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, [0.3, 0.1, 0.2])
    assert exc.value.args[0] == attribute + ' should be monotonically increasing'


def test_range_nu_valid1():
    o = OpticalProperties()
    o.nu = [0.1, 0.5, 0.8]


def test_range_nu_invalid1():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.nu = [0., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_nu_invalid2():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.nu = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_mu_valid1():
    o = OpticalProperties()
    o.mu = [-0.5, 0., 0.5]


def test_range_mu_valid2():
    o = OpticalProperties()
    o.mu = [-1., 0., 1.]


def test_range_mu_invalid1():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.mu = [-1.3, 0., 1.]
    assert exc.value.args[0] == 'mu should be in the range [-1:1]'


def test_range_mu_invalid2():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.mu = [-1., 0., 1.3]
    assert exc.value.args[0] == 'mu should be in the range [-1:1]'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_list(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    setattr(o, attribute, [[1., 2.], [0., 1.], [3., 4.]])


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_array(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    setattr(o, attribute, np.ones((3, 2)))


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_type1(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, 'hello')
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_type2(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, 2.123)
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_shape1(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, [1., 2., 3.])
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_shape2(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    o.mu = [-0.5, 0.5]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, np.ones((4, 5)))
    assert exc.value.args[0] == attribute + ' has an incorrect shape: (4, 5), but expected (3, 2)'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_order1(attribute):
    o = OpticalProperties()
    o.nu = [0.1, 0.2, 0.3]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, np.ones((3, 2)))
    assert exc.value.args[0] == 'mu needs to be set before ' + attribute


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_order2(attribute):
    o = OpticalProperties()
    o.mu = [-0.5, 0.5]
    with pytest.raises(ValueError) as exc:
        setattr(o, attribute, np.ones((3, 2)))
    assert exc.value.args[0] == 'nu needs to be set before ' + attribute
