from __future__ import print_function, division

import matplotlib.pyplot as plt

import numpy as np
from numpy.testing import assert_allclose, assert_array_almost_equal_nulp

from astropy.tests.helper import pytest

from ..optical_properties import OpticalProperties
from ...util.constants import c
from ...util.functions import virtual_file


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


@pytest.mark.parametrize(('attribute'), ['nu', 'mu'])
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


def test_range_chi_valid1():
    o = OpticalProperties()
    o.chi = [0.1, 0.5, 0.8]


def test_range_chi_valid2():
    o = OpticalProperties()
    o.chi = [0., 0.5, 0.8]


def test_range_chi_invalid1():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.chi = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'chi should be positive'


def test_range_albedo_valid1():
    o = OpticalProperties()
    o.albedo = [0., 0.5, 1.]


def test_range_albedo_invalid1():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.albedo = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'albedo should be in the range [0:1]'


def test_range_albedo_invalid2():
    o = OpticalProperties()
    with pytest.raises(ValueError) as exc:
        o.albedo = [0., 0.5, 1.1]
    assert exc.value.args[0] == 'albedo should be in the range [0:1]'


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
    assert exc.value.args[0] == attribute + ' has an incorrect shape: (4, 5) but expected (3, 2)'


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


def test_extrapolate_inner_range():
    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.extrapolate_nu(1e9, 2e9)
    assert o.nu[0] == 1.e8 and o.nu[-1] == 1.e10


def test_extrapolate_upper():
    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.extrapolate_nu(1e9, 1e11)
    assert o.nu[0] == 1.e8 and o.nu[-1] == 1.e11


def test_extrapolate_lower():
    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.extrapolate_nu(1e7, 1e9)
    assert o.nu[0] == 1.e7 and o.nu[-1] == 1.e10


def test_extrapolate_both():
    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.extrapolate_nu(1e7, 1e11)
    assert o.nu[0] == 1.e7 and o.nu[-1] == 1.e11


def test_extrapolate_wav():
    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.extrapolate_wav(1., 1.e20)
    assert_array_almost_equal_nulp(o.nu[0], c / 1.e16, 2)
    assert_array_almost_equal_nulp(o.nu[-1], c / 1.e-4, 2)



def test_io():

    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()

    f = virtual_file()
    o.to_hdf5_group(f)
    o_new = OpticalProperties()
    o_new.from_hdf5_group(f)

    assert_allclose(o.nu, o_new.nu)
    assert_allclose(o.chi, o_new.chi)
    assert_allclose(o.albedo, o_new.albedo)
    assert_allclose(o.mu, o_new.mu)
    assert_allclose(o.P1, o_new.P1)
    assert_allclose(o.P2, o_new.P2)
    assert_allclose(o.P3, o_new.P3)
    assert_allclose(o.P4, o_new.P4)
    assert o.hash() == o_new.hash()


def test_plot():

    # Just check that plot runs without crashing

    fig = plt.figure()

    o = OpticalProperties()
    o.nu = np.logspace(8., 10., 100)
    o.albedo = np.repeat(0.5, 100)
    o.chi = np.ones(100)
    o.mu = [-1., 1.]
    o.initialize_scattering_matrix()
    o.plot(fig, [321,322,323,324,325,326])

    plt.close(fig)
