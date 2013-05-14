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


# TODO: implement tests for RunConf and ImageConf (and subclasses)
