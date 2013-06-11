from astropy.tests.helper import pytest
import numpy as np

from ..conf_files import PeeledImageConf
from ...util.functions import virtual_file


# Viewing angles


def test_viewing_angles_tuple():
    i = PeeledImageConf()
    i.set_viewing_angles([1., 2., 3.], [4., 5., 6.])
    i._write_viewing_angles(virtual_file())


def test_viewing_angles_list():
    i = PeeledImageConf()
    i.set_viewing_angles([1., 2., 3.], [4., 5., 6.])
    i._write_viewing_angles(virtual_file())


def test_viewing_angles_array():
    i = PeeledImageConf()
    i.set_viewing_angles(np.array([1., 2., 3.]), np.array([4., 5., 6.]))
    i._write_viewing_angles(virtual_file())


def test_viewing_angles_theta_type_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_viewing_angles('test', [4., 5.])
    assert exc.value.args[0] == 'theta should be a 1-D sequence'


def test_viewing_angles_phi_type_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_viewing_angles([4., 5.], 'test')
    assert exc.value.args[0] == 'phi should be a 1-D sequence'


def test_viewing_angles_theta_dims_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_viewing_angles([[1., 2., 3.], [4., 5., 6.]], [4., 5.])
    assert exc.value.args[0] == 'theta should be a 1-D sequence'


def test_viewing_angles_phi_dims_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_viewing_angles([4., 5.], [[1., 2., 3.], [4., 5., 6.]])
    assert exc.value.args[0] == 'phi should be a 1-D sequence'


def test_viewing_angles_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_viewing_angles([1., 2., 3.], [4., 5.])
    assert exc.value.args[0] == 'Length of theta and phi arrays do not match'


# Inside observer position


@pytest.mark.parametrize(('position'), [None,
                                        (1., 2., 3.),
                                        [1., 2., 3.],
                                        np.array([1., 2., 3.])])
def test_inside_observer(position):
    i = PeeledImageConf()
    i.set_inside_observer(position)
    if position is not None:
        i._write_inside_observer(virtual_file())


def test_inside_observer_phi_type_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_inside_observer('test')
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


def test_inside_observer_phi_dims_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_inside_observer([[1., 2., 3.], [4., 5., 6.]])
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


def test_inside_observer_phi_len_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_inside_observer([1., 2., 3., 4.])
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


# Peeloff origin


@pytest.mark.parametrize(('position'), [None, (1., 2., 3.), [1., 2., 3.], np.array([1., 2., 3.])])
def test_peeloff_origin(position):
    i = PeeledImageConf()
    i.set_peeloff_origin(position)
    if position is not None:
        i._write_peeloff_origin(virtual_file())


def test_peeloff_origin_phi_type_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_peeloff_origin('test')
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


def test_peeloff_origin_phi_dims_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_peeloff_origin([[1., 2., 3.], [4., 5., 6.]])
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


def test_peeloff_origin_phi_len_mismatch():
    i = PeeledImageConf()
    with pytest.raises(ValueError) as exc:
        i.set_peeloff_origin([1., 2., 3., 4.])
    assert exc.value.args[0] == 'position should be a 1-D sequence with 3 elements'


# Main write method


def test_main_missing():
    i = PeeledImageConf()
    with pytest.raises(Exception) as exc:
        i._write_viewing_info(virtual_file())
    assert exc.value.args[0] == 'Need to specify either observer position, or viewing angles'


def test_main_inside():
    i = PeeledImageConf()
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(10., -10., -5., 5.)
    i._write_viewing_info(virtual_file())


def test_main_invalid_longitude():
    i = PeeledImageConf()
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(-10., 10., -5., 5.)
    with pytest.raises(Exception) as exc:
        i._write_viewing_info(virtual_file())
    assert exc.value.args[0] == 'longitudes should increase towards the left for inside observers'


def test_main_invalid_depth1():
    i = PeeledImageConf()
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(10., -10., -5., 5.)
    i.set_depth(-1., 1.)
    with pytest.raises(ValueError) as exc:
        i._write_viewing_info(virtual_file())
    assert exc.value.args[0] == 'Lower limit of depth should be positive for inside observer'


def test_main_invalid_depth2():
    i = PeeledImageConf()
    i.set_inside_observer((0., 0., 0.))
    i.set_image_limits(10., -10., -5., 5.)
    i.set_depth(0., -1.)
    with pytest.raises(ValueError) as exc:
        i._write_viewing_info(virtual_file())
    assert exc.value.args[0] == 'Upper limit of depth should be positive for inside observer'


def test_main_invalid_inside_and_peeloff():
    i = PeeledImageConf()
    i.set_inside_observer((0., 0., 0.))
    i.set_peeloff_origin((0., 0., 0.))
    with pytest.raises(Exception) as exc:
        i._write_viewing_info(virtual_file())
    assert exc.value.args[0] == 'Cannot specify inside observer and peeloff origin at the same time'
