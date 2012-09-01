import pytest
import numpy as np

from ..surface_scattering_properties import SurfaceScatteringProperties


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


def test_init():
    SurfaceScatteringProperties()

BASE_ATTRIBUTES = ['nu', 'mu0', 'mu', 'psi']

# TYPES


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_list(attribute):
    p = SurfaceScatteringProperties()
    setattr(p, attribute, [0.1, 0.2, 0.3])


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_array(attribute):
    p = SurfaceScatteringProperties()
    setattr(p, attribute, np.array([0.1, 0.2, 0.3]))


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_invalid_type1(attribute):
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        setattr(p, attribute, 'hello')
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_invalid_type2(attribute):
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        setattr(p, attribute, 0.5)
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_invalid_shape1(attribute):
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        setattr(p, attribute, [[0., 1.], [0.5, 1.]])
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_invalid_shape2(attribute):
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        setattr(p, attribute, np.array([[0., 1.], [0.5, 1.]]))
    assert exc.value.args[0] == attribute + ' should be a 1-D sequence'

# TYPES (ALBEDO)


def test_set_albedo_list():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.albedo = [0.1, 0.2, 0.3]


def test_set_albedo_array():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.albedo = np.array([0.1, 0.2, 0.3])


def test_set_albedo_invalid_type1():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    with pytest.raises(ValueError) as exc:
        p.albedo = 'hello'
    assert exc.value.args[0] == 'albedo should be a 1-D sequence'


def test_set_albedo_invalid_type2():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    with pytest.raises(ValueError) as exc:
        p.albedo = 0.5
    assert exc.value.args[0] == 'albedo should be a 1-D sequence'


def test_set_albedo_invalid_shape1():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    with pytest.raises(ValueError) as exc:
        p.albedo = [[0., 1.], [0.5, 1.]]
    assert exc.value.args[0] == 'albedo should be a 1-D sequence'


def test_set_albedo_invalid_shape2():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    with pytest.raises(ValueError) as exc:
        p.albedo = np.array([[0., 1.], [0.5, 1.]])
    assert exc.value.args[0] == 'albedo should be a 1-D sequence'

# ORDER


@pytest.mark.parametrize(('attribute'), BASE_ATTRIBUTES)
def test_set_vector_invalid_order(attribute):
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        setattr(p, attribute, [0.3, 0.1, 0.2])
    assert exc.value.args[0] == attribute + ' should be monotonically increasing'

# RANGES


def test_range_nu_valid1():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]


def test_range_nu_invalid1():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.nu = [0., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_nu_invalid2():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.nu = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'nu should be strictly positive'


def test_range_albedo_valid1():
    p = SurfaceScatteringProperties()
    p.nu = [1., 2., 3.]
    p.albedo = [0., 0.5, 1.]


def test_range_albedo_nonu():
    p = SurfaceScatteringProperties()
    with pytest.raises(Exception) as exc:
        p.albedo = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'nu has to be defined first'


def test_range_albedo_invalid1():
    p = SurfaceScatteringProperties()
    p.nu = [1., 2., 3.]
    with pytest.raises(ValueError) as exc:
        p.albedo = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'albedo should be in the range [0:1]'


def test_range_albedo_invalid2():
    p = SurfaceScatteringProperties()
    p.nu = [1., 2., 3.]
    with pytest.raises(ValueError) as exc:
        p.albedo = [0., 0.5, 1.1]
    assert exc.value.args[0] == 'albedo should be in the range [0:1]'


def test_range_mu0_valid1():
    p = SurfaceScatteringProperties()
    p.mu0 = [0., 0.5, 1.]


def test_range_mu0_invalid1():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.mu0 = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'mu0 should be in the range [0:1]'


def test_range_mu0_invalid2():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.mu0 = [0., 0.5, 1.1]
    assert exc.value.args[0] == 'mu0 should be in the range [0:1]'


def test_range_mu_valid1():
    p = SurfaceScatteringProperties()
    p.mu = [0., 0.5, 1.]


def test_range_mu_invalid1():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.mu = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'mu should be in the range [0:1]'


def test_range_mu_invalid2():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.mu = [0., 0.5, 1.1]
    assert exc.value.args[0] == 'mu should be in the range [0:1]'


def test_range_psi_valid1():
    p = SurfaceScatteringProperties()
    p.psi = [0., 0.5, 1.]


def test_range_psi_invalid1():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.psi = [-1., 0.5, 0.8]
    assert exc.value.args[0] == 'psi should be in the range [0:2pi]'


def test_range_psi_invalid2():
    p = SurfaceScatteringProperties()
    with pytest.raises(ValueError) as exc:
        p.psi = [0., 0.5, 10.]
    assert exc.value.args[0] == 'psi should be in the range [0:2pi]'


# BRDF

def test_set_brdf():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu0 = [0., 1.]
    p.mu = [0., 0.33, 0.66, 1.0]
    p.psi = [0., 1., 2., 3., 4.]
    p.brdf = np.zeros((3, 2, 4, 5))


def test_set_brdf_nonu():
    p = SurfaceScatteringProperties()
    p.mu0 = [0., 1.]
    p.mu = [0., 0.33, 0.66, 1.0]
    p.psi = [0., 1., 2., 3., 4.]
    with pytest.raises(Exception) as exc:
        p.brdf = np.zeros((3, 2, 4, 5))
    assert exc.value.args[0] == 'nu has to be defined first'


def test_set_brdf_nomu0():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu = [0., 0.33, 0.66, 1.0]
    p.psi = [0., 1., 2., 3., 4.]
    with pytest.raises(Exception) as exc:
        p.brdf = np.zeros((3, 2, 4, 5))
    assert exc.value.args[0] == 'mu0 has to be defined first'


def test_set_brdf_nomu():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu0 = [0., 1.]
    p.psi = [0., 1., 2., 3., 4.]
    with pytest.raises(Exception) as exc:
        p.brdf = np.zeros((3, 2, 4, 5))
    assert exc.value.args[0] == 'mu has to be defined first'


def test_set_brdf_nopsi():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu0 = [0., 1.]
    p.mu = [0., 0.33, 0.66, 1.0]
    with pytest.raises(Exception) as exc:
        p.brdf = np.zeros((3, 2, 4, 5))
    assert exc.value.args[0] == 'psi has to be defined first'


@pytest.mark.parametrize('value', ['hello', np.zeros((2,)), np.zeros((4, 2)),
                         np.zeros((3, 2, 3)), [[[[1.]]]]])
def test_set_brdf_invalid_type(value):
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu0 = [0., 1.]
    p.mu = [0., 0.33, 0.66, 1.0]
    p.psi = [0., 1., 2., 3., 4.]
    with pytest.raises(ValueError) as exc:
        p.brdf = value
    assert exc.value.args[0] == 'brdf should be a 4-d Numpy array'


def test_set_brdf_invalid_shape():
    p = SurfaceScatteringProperties()
    p.nu = [0.1, 0.5, 0.8]
    p.mu0 = [0., 1.]
    p.mu = [0., 0.33, 0.66, 1.0]
    p.psi = [0., 1., 2., 3., 4.]
    with pytest.raises(ValueError) as exc:
        p.brdf = np.zeros((2, 4, 3, 3))
    assert exc.value.args[0] == 'brdf has an incorrect shape: (2, 4, 3, 3) but expected (3, 2, 4, 5)'


def test_io_roundtrip(tmpdir):

    p1 = SurfaceScatteringProperties()
    p1.nu = [0.1, 0.5, 0.8]
    p1.mu0 = [0., 1.]
    p1.mu = [0., 0.33, 0.66, 1.0]
    p1.psi = [0., 1., 2., 3., 4.]
    p1.albedo = [0.2, 0.3, 0.4]
    p1.brdf = np.arange(120).reshape(3, 2, 4, 5)

    filename = str(tmpdir.join('test_round_trip'))
    p1.write(filename)
    p2 = SurfaceScatteringProperties(filename)

    assert np.all(p1.nu == p2.nu)
    assert np.all(p1.mu0 == p2.mu0)
    assert np.all(p1.mu == p2.mu)
    assert np.all(p1.psi == p2.psi)
    assert np.all(p1.albedo == p2.albedo)
    assert np.all(p1.brdf == p2.brdf)
