from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np

from ..optical_properties import OpticalProperties
from ..emissivities import Emissivities
from ...util.functions import virtual_file


def test_init():
    Emissivities()

VECTOR_ATTRIBUTES = ['nu', 'var']
ARRAY_ATTRIBUTES = ['jnu']


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_list(attribute):
    e = Emissivities()
    setattr(e, attribute, [0.1, 0.2, 0.3])


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_array(attribute):
    e = Emissivities()
    setattr(e, attribute, np.array([0.1, 0.2, 0.3]))


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_type1(attribute):
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, 'hello')
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_type2(attribute):
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, 0.5)
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_shape1(attribute):
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, [[0., 1.], [0.5, 1.]])
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_shape2(attribute):
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, np.array([[0., 1.], [0.5, 1.]]))
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), VECTOR_ATTRIBUTES)
def test_set_vector_invalid_order(attribute):
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, [0.3, 0.1, 0.2])
    assert exc.value.args[0] == attribute + ' should be monotonically increasing'


def test_range_nu_valid1():
    e = Emissivities()
    e.nu = [0.1, 0.5, 0.8]


def test_range_nu_invalid1():
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        e.nu = [0., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_nu_invalid2():
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        e.nu = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_var_valid1():
    e = Emissivities()
    e.var = [0.1, 0.5, 0.8]


def test_range_var_valid2():
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        e.var = [0., 0.5, 0.8]
    assert exc.value.args[0] == 'var should be strictly positive'


def test_range_var_invalid1():
    e = Emissivities()
    with pytest.raises(ValueError) as exc:
        e.var = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'var should be strictly positive'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_list(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    setattr(e, attribute, [[1., 2.], [0., 1.], [3., 4.]])


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_array(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    setattr(e, attribute, np.ones((3, 2)))


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_type1(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, 'hello')
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_type2(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, 2.123)
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_shape1(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, [1., 2., 3.])
    assert exc.value.args[0] == attribute + ' should be a 2-D array'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_shape2(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, np.ones((4, 5)))
    assert exc.value.args[0] == attribute + ' has an incorrect shape: (4, 5) but expected (3, 2)'


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_order1(attribute):
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, np.ones((3, 2)))
    assert exc.value.args[0] == 'var needs to be set before ' + attribute


@pytest.mark.parametrize(('attribute'), ARRAY_ATTRIBUTES)
def test_set_array_invalid_order2(attribute):
    e = Emissivities()
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        setattr(e, attribute, np.ones((3, 2)))
    assert exc.value.args[0] == 'nu needs to be set before ' + attribute


def test_set_jnu_valid1():
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    e.jnu = np.ones((3, 2))


def test_set_jnu_valid2():
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    e.jnu = np.zeros((3, 2))


def test_set_jnu_invalid1():
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    with pytest.raises(ValueError) as exc:
        e.jnu = np.ones((3, 2)) * -1.
    assert exc.value.args[0] == 'jnu should be positive'


def test_lte_emissivities():
    o = OpticalProperties()
    o.nu = np.array([1.e8, 1.e16])
    o.chi = np.array([1.e-2, 1])
    o.albedo = np.array([0., 0.5])
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    e.set_lte(o, n_temp=10, temp_min=1., temp_max=1000.)


def test_io():
    o = OpticalProperties()
    o.nu = np.array([1.e8, 1.e16])
    o.chi = np.array([1.e-2, 1])
    o.albedo = np.array([0., 0.5])
    e = Emissivities()
    e.nu = [0.1, 0.2, 0.3]
    e.var = [0.1, 1.]
    e.set_lte(o, n_temp=10, temp_min=1., temp_max=1000.)
    f = virtual_file()
    e.to_hdf5_group(f)
    e_new = Emissivities()
    e_new.from_hdf5_group(f)
    assert e.var_name == e_new.var_name
    assert np.all(e.nu == e_new.nu)
    assert np.all(e.var == e_new.var)
    assert np.all(e.jnu == e_new.jnu)
