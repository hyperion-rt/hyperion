from __future__ import print_function, division

import string
import random

import h5py
from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import Model


def random_id(length=32):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


@pytest.mark.parametrize(('value'), [True, False])
def test_io_monochromatic(value):
    m1 = Model()
    if value:
        m1.set_monochromatic(value, frequencies=[1, 4.5, 7.7])
    else:
        m1.set_monochromatic(value)
    v = virtual_file()
    m1._write_monochromatic(v)
    m2 = Model()
    m2._read_monochromatic(v)
    assert m1._monochromatic == m2._monochromatic
    if value:
        assert_equal(m1._frequencies, m2._frequencies)
