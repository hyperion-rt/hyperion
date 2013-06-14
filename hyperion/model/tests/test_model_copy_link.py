from __future__ import print_function, division

import os
import shutil
import tempfile

import h5py
from astropy.tests.helper import pytest
import numpy as np

from ...util.functions import random_id
from .test_helpers import get_test_model_noimaging, get_test_dust


class TestWriteDustCopy(object):

    def setup_method(self, method):

        self.dust = get_test_dust()

        self.density = np.array([[[1.]]])

        self.model = get_test_model_noimaging()

    def test_copy_filename(self, tmpdir):
        dust_file = tmpdir.join(random_id()).strpath
        self.dust.write(dust_file)
        self.model.add_density_grid(self.density, dust_file)
        self.model.write(tmpdir.join(random_id()).strpath, copy=True)
        self.model.run(tmpdir.join(random_id()).strpath)

    def test_copy_object(self, tmpdir):
        self.model.add_density_grid(self.density, self.dust)
        self.model.write(tmpdir.join(random_id()).strpath, copy=True)
        self.model.run(tmpdir.join(random_id()).strpath)

    def test_link_filename(self, tmpdir):
        dust_file = tmpdir.join(random_id()).strpath
        self.dust.write(dust_file)
        self.model.add_density_grid(self.density, dust_file)
        self.model.write(tmpdir.join(random_id()).strpath, copy=False)
        self.model.run(tmpdir.join(random_id()).strpath)

    def test_link_object(self, tmpdir):
        self.model.add_density_grid(self.density, self.dust)
        with pytest.raises(ValueError) as e:
            self.model.write(tmpdir.join(random_id()).strpath, copy=False)
        assert e.value.args[0] == 'Dust properties are not located in a file, so cannot link. Use copy=True or write the dust properties to a file first'


@pytest.mark.parametrize(('write_copy'), [True, False])
def test_input_link(write_copy, tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.set_copy_input(False)
    model.write(input_file, copy=write_copy)

    # Check that copy parameter is there
    f = h5py.File(input_file, 'r')
    assert f.attrs['copy_input'].decode('utf-8') == 'no'
    f.close()

    # Run the model
    model.run(output_file)

    # Check that attributes from input can still be accessed, then check that
    # 'Input' is a link
    f = h5py.File(output_file, 'r')
    assert f['Input'].attrs['copy_input'].decode('utf-8') == 'no'
    assert f.file != f['Input'].file
    f.close()


@pytest.mark.parametrize(('write_copy'), [True, False])
def test_input_copy(write_copy, tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    model = get_test_model_noimaging()
    model.set_copy_input(True)
    model.write(input_file, copy=write_copy)

    # Check that copy parameter is there
    f = h5py.File(input_file, 'r')
    assert f.attrs['copy_input'].decode('utf-8') == 'yes'
    f.close()

    # Run the model
    model.run(output_file)

    # Check that attributes from input can still be accessed, then check that
    # 'Input' is a not a link
    f = h5py.File(output_file, 'r')
    assert f['Input'].attrs['copy_input'].decode('utf-8') == 'yes'
    assert f.file == f['Input'].file
    f.close()
