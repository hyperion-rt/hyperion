from __future__ import print_function, division

from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import Model
from ...util.functions import virtual_file
from .test_helpers import random_id


def test_io_minimal(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = Model()
    m1.set_cartesian_grid([-1., 1.],[-1., 1.],[-1., 1.])
    m1.set_n_initial_iterations(0)
    m1.set_n_photons(imaging=10)
    m1.write(filename)

    m2 = Model.read(filename)


def test_io_minimal_mono(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = Model()
    m1.set_cartesian_grid([-1., 1.],[-1., 1.],[-1., 1.])
    m1.set_n_initial_iterations(0)
    m1.set_monochromatic(True, frequencies=[1, 4.5, 7.7])
    m1.set_n_photons(imaging_sources=10, imaging_dust=10)
    m1.write(filename)

    m2 = Model.read(filename)

    assert m1._monochromatic == m2._monochromatic
    assert_equal(m1._frequencies, m2._frequencies)


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
