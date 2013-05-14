# Test that parameters are preserved when written out and read in again

from __future__ import print_function, division

import string
import random

import h5py
import numpy as np
from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import OutputConf, RunConf, ImageConf, BinnedImageConf, PeeledImageConf


def random_id(length=32):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


def test_io_output_conf():
    o1 = OutputConf()
    o1.output_density = 'last'
    o1.output_density_diff = 'none'
    o1.output_specific_energy = 'all'
    o1.output_n_photons = 'last'
    v = virtual_file()
    o1.write(v)
    o2 = OutputConf.read(v)
    assert o2.output_density == o1.output_density
    assert o2.output_density_diff == o1.output_density_diff
    assert o2.output_specific_energy == o1.output_specific_energy
    assert o2.output_n_photons == o1.output_n_photons


# TODO: implement tests for RunConf and ImageConf (and subclasses)
