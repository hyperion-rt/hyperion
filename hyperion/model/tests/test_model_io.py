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


def test_io_monochromatic_full(tmpdir):
    """
    Regression for issues described in #77 - namely that when reading and
    writing out a model with monochromatic settings, images/SEDs would not work
    properly (raising an exception).
    """

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    i = m.add_peeled_images(image=True)
    i.set_viewing_angles([45.], [45.])
    i.set_image_size(256, 256)
    i.set_image_limits(-1., 1., -1., 1.)

    m.set_monochromatic(True, wavelengths=[0.1, 1., 10.])

    m.set_n_photons(initial=1000, imaging_sources=1000, imaging_dust=1000)

    m.write(filename1)

    m2 = Model.read(filename1)
    m2.write(filename2)

def test_io_monochromatic_full_2(tmpdir):
    """
    Regression for issues described in #77 - namely that when reading and
    writing out a model with monochromatic settings, images/SEDs would not work
    properly (raising an exception).

    Changing the order to setting before the monochromatic mode is set changes
    how things are handled internally, hence the similar test.
    """

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    m.set_monochromatic(True, wavelengths=[0.1, 1., 10.])

    i = m.add_peeled_images(image=True)
    i.set_viewing_angles([45.], [45.])
    i.set_image_size(256, 256)
    i.set_image_limits(-1., 1., -1., 1.)

    m.set_n_photons(initial=1000, imaging_sources=1000, imaging_dust=1000)

    m.write(filename1)

    m2 = Model.read(filename1)
    m2.write(filename2)


def test_io_voronoi(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    x = [1., 3., 2., 5., 2.]
    y = [3., 2., 8., 3., 6.]
    z = [4., 3., 1., 1., 2.]

    m1 = Model()
    m1.set_voronoi_grid(x, y, z)
    m1.set_n_initial_iterations(0)
    m1.set_n_photons(imaging=10)
    m1.write(filename)

    m2 = Model.read(filename)
