import os
import shutil
import tempfile

import numpy as np

from astropy.tests.helper import pytest
from numpy.testing import assert_array_almost_equal_nulp

from .test_helpers import get_test_model_noimaging, random_id, get_test_dust


def setup_module(module):
    module.tmpdir = tempfile.mkdtemp()
    module.dust_file = os.path.join(module.tmpdir, random_id())
    dust = get_test_dust()
    dust.write(module.dust_file)
    module.density = np.array([[[1.]]])


def teardown_module(module):
    shutil.rmtree(module.tmpdir)


def test_minimum_temperature_scalar(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature(10.)
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('temperature')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 10., 10)


def test_minimum_temperature_scalar_list(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature([10.])
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('temperature')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 10., 10)


def test_minimum_temperature_scalar_invalid1():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_temperature(-10.)
    assert exc.value.args[0] == 'temperature should be positive'


def test_minimum_temperature_scalar_invalid2():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_temperature('a')
    assert exc.value.args[0] == 'temperature should be a numerical value'


def test_minimum_temperature_scalar_invalid3():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_temperature([-10.])
    assert exc.value.args[0] == 'temperature should be positive'


def test_minimum_temperature_scalar_invalid4():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_temperature(['a'])
    assert exc.value.args[0] == 'temperature should be a numerical value'


def test_minimum_temperature_scalar_invalid5(tmpdir):

    input_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature([10., 10.])
    with pytest.raises(Exception) as exc:
        model.write(input_file)
    assert exc.value.args[0] == 'Number of minimum_temperature values should match number of dust types'


def test_minimum_temperature_scalar_2(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature(10.)
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('temperature')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 10., 10)
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 10., 10)


def test_minimum_temperature_scalar_list_2(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature([10., 8.])
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('temperature')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 10., 10)
    assert_array_almost_equal_nulp(t[1][0, 0, 0], 8., 10)


def test_minimum_temperature_scalar_list_2_invalid(tmpdir):

    input_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_temperature([10., 8., 6.])
    with pytest.raises(Exception) as exc:
        model.write(input_file)
    assert exc.value.args[0] == 'Number of minimum_temperature values should match number of dust types'


def test_minimum_specific_energy_scalar(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy(2.)
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('specific_energy')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 2., 10)


def test_minimum_specific_energy_scalar_list(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy([2.])
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('specific_energy')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 2., 10)


def test_minimum_specific_energy_scalar_invalid1():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_specific_energy(-10.)
    assert exc.value.args[0] == 'specific_energy should be positive'


def test_minimum_specific_energy_scalar_invalid2():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_specific_energy('a')
    assert exc.value.args[0] == 'specific_energy should be a numerical value'


def test_minimum_specific_energy_scalar_invalid3():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_specific_energy([-10.])
    assert exc.value.args[0] == 'specific_energy should be positive'


def test_minimum_specific_energy_scalar_invalid4():

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    with pytest.raises(Exception) as exc:
        model.set_minimum_specific_energy(['a'])
    assert exc.value.args[0] == 'specific_energy should be a numerical value'


def test_minimum_specific_energy_scalar_invalid5(tmpdir):

    input_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy([2., 2.])
    with pytest.raises(Exception) as exc:
        model.write(input_file)
    assert exc.value.args[0] == 'Number of minimum_specific_energy values should match number of dust types'


def test_minimum_specific_energy_scalar_2(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy(2.)
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('specific_energy')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 2., 10)
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 2., 10)


def test_minimum_specific_energy_scalar_list_2(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy([2., 3.])
    model.write(input_file)
    out = model.run(output_file)
    t = out.get_physical_grid('specific_energy')
    assert_array_almost_equal_nulp(t[0][0, 0, 0], 2., 10)
    assert_array_almost_equal_nulp(t[1][0, 0, 0], 3., 10)


def test_minimum_specific_energy_scalar_list_2_invalid(tmpdir):

    input_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.add_density_grid(density, dust_file)
    model.add_density_grid(density, dust_file)
    model.set_minimum_specific_energy([2., 3., 4.])
    with pytest.raises(Exception) as exc:
        model.write(input_file)
    assert exc.value.args[0] == 'Number of minimum_specific_energy values should match number of dust types'
