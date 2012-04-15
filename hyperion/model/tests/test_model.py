from __future__ import print_function, division

from .. import Model
from .test_helpers import random_filename


def test_basic():

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())
