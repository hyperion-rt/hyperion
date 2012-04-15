from __future__ import print_function, division

import os
import tempfile
import numpy as np

from ...dust import IsotropicSphericalDust
from .. import Model
from ...util.functions import random_id


def random_filename():
    return os.path.join(tempfile.mkdtemp(), random_id())


def get_test_dust():
    dust = IsotropicSphericalDust([1.e-2, 1.e5], [1., 1.], [0.5, 0.5])
    dust.emissivities.set_lte(dust.optical_properties, n_temp=10, temp_min=0.1, temp_max=1600.)
    return dust


def get_test_model_noimaging():

    model = Model()
    model.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    model.set_n_photons(initial=1, imaging=0)
    model.set_n_initial_iterations(1)

    source = model.add_point_source()
    source.luminosity = 1.
    source.temperature = 1000.

    return model
