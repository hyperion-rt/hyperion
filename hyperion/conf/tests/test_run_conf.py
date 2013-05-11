from astropy.tests.helper import pytest
from ..conf_files import RunConf


@pytest.mark.parametrize(('value'), [0, 0., 0.2, 0.523, 1., 1])
def test_propagation_check_frequency(value):
    r = RunConf()
    r.set_propagation_check_frequency(value)


@pytest.mark.parametrize(('value'), [[1, 2, 3], 'hello', (1, 2)])
def test_propagation_check_frequency_invalid1(value):
    r = RunConf()
    with pytest.raises(TypeError) as exc:
        r.set_propagation_check_frequency(value)
    assert exc.value.args[0] == "frequency should be a scalar value"


@pytest.mark.parametrize(('value'), [-1., -0.3, 1.3, 1e20])
def test_propagation_check_frequency_invalid2(value):
    r = RunConf()
    with pytest.raises(ValueError) as exc:
        r.set_propagation_check_frequency(value)
    assert exc.value.args[0] == "frequency should be between 0 and 1"
