import pytest
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


def test_set_convergence_flag():
    r = RunConf()
    r.set_convergence(True, percentile=99., absolute=2., relative=1.02)
    assert r.check_convergence is True
    r = RunConf()
    r.set_convergence(False)
    assert r.check_convergence is False


@pytest.mark.parametrize(('value'), [1e-3, 1, 50., 100, 1e300, float('inf')])
def test_max_peeloff_optical_depth(value):
    r = RunConf()
    r.set_max_peeloff_optical_depth(value)
    assert r.peeloff_tau_max == value


@pytest.mark.parametrize(('value'), [[1, 2, 3], 'hello', (1, 2)])
def test_max_peeloff_optical_depth_invalid1(value):
    r = RunConf()
    with pytest.raises(TypeError) as exc:
        r.set_max_peeloff_optical_depth(value)
    assert exc.value.args[0] == "tau_max should be a scalar value"


@pytest.mark.parametrize(('value'), [0., -1., float('nan')])
def test_max_peeloff_optical_depth_invalid2(value):
    r = RunConf()
    with pytest.raises(ValueError) as exc:
        r.set_max_peeloff_optical_depth(value)
    assert exc.value.args[0] == "tau_max should be strictly positive"
