# Test that parameters are preserved when written out and read in again

from __future__ import print_function, division

import string
import random
from itertools import product

import h5py
import numpy as np
from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import OutputConf, RunConf, ImageConf, BinnedImageConf, PeeledImageConf


def random_id(length=32):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


@pytest.mark.parametrize(('attribute', 'value'),
                         list(product(['output_density', 'output_density_diff',
                         'output_specific_energy', 'output_n_photons'],
                         ['none', 'last', 'all'])))
def test_io_output_conf(attribute, value):
    o1 = OutputConf()
    setattr(o1, attribute, value)
    v = virtual_file()
    o1.write(v)
    o2 = OutputConf.read(v)
    assert getattr(o2, attribute) == value


def test_io_image_conf():
    i1 = ImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_output_bytes(4)
    i1.set_track_origin('basic')
    i1.set_uncertainties(True)
    v = virtual_file()
    i1.write(v)
    i2 = ImageConf.read(v)
    assert i2.n_x == i1.n_x
    assert i2.n_y == i1.n_y
    assert i2.xmin == i1.xmin
    assert i2.xmax == i1.xmax
    assert i2.ymin == i1.ymin
    assert i2.ymax == i1.ymax
    assert i2.n_ap == i1.n_ap
    assert i2.ap_min == i1.ap_min
    assert i2.ap_max == i1.ap_max
    assert i2.n_wav == i1.n_wav
    assert i2.wav_min == i1.wav_min
    assert i2.wav_max == i1.wav_max
    assert i2.io_bytes == i1.io_bytes
    assert i2.track_origin == i1.track_origin
    assert i2.uncertainties == i1.uncertainties


def test_io_binned_image_conf():
    i1 = BinnedImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_bins(76, 22)
    v = virtual_file()
    i1.write(v)
    i2 = BinnedImageConf.read(v)
    assert i2.n_theta == i1.n_theta
    assert i2.n_phi == i1.n_phi


def test_io_peeled_image_conf():
    i1 = PeeledImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_angles([1.,2.,3], [4.,5.,6.])
    i1.set_peeloff_origin([2.2, 3.3, 7.6])
    i1.set_ignore_optical_depth(True)
    i1.set_depth(-1.7, 6.2)
    v = virtual_file()
    i1.write(v)
    i2 = PeeledImageConf.read(v)
    for i in range(len(i2.viewing_angles)):
        assert i2.viewing_angles[i][0] == i1.viewing_angles[i][0]
        assert i2.viewing_angles[i][1] == i1.viewing_angles[i][1]
    assert_equal(i2.peeloff_origin, i1.peeloff_origin)
    assert i2.ignore_optical_depth == i1.ignore_optical_depth
    assert i2.d_min == i1.d_min
    assert i2.d_max == i1.d_max


def test_io_peeled_image_conf_inside():
    i1 = PeeledImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, -4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_angles([1.,2.,3], [4.,5.,6.])
    i1.set_inside_observer([7.,8.,9.])
    i1.set_ignore_optical_depth(True)
    i1.set_depth(1.7, 6.2)
    v = virtual_file()
    i1.write(v)
    i2 = PeeledImageConf.read(v)
    for i in range(len(i2.viewing_angles)):
        print(i2.viewing_angles[i][0])
        assert i2.viewing_angles[i][0] == i1.viewing_angles[i][0]
        assert i2.viewing_angles[i][1] == i1.viewing_angles[i][1]
    assert_equal(i2.inside_observer, i1.inside_observer)
    assert i2.ignore_optical_depth == i1.ignore_optical_depth
    assert i2.d_min == i1.d_min
    assert i2.d_max == i1.d_max

# TODO: implement tests for RunConf and ImageConf (and subclasses)
