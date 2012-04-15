from __future__ import print_function, division

import pytest

from .. import AnalyticalYSOModel
from ...util.constants import msun, rsun, lsun, tsun, au, pc
from ...util.functions import random_filename
from .test_helpers import get_test_dust


def test_full_analytical_yso():

    dust = get_test_dust()

    m = AnalyticalYSOModel()

    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun

    d = m.add_flared_disk()
    d.mass = 1.e-2 * msun
    d.rmax = 300. * au
    d.rmin = 0.1 * au
    d.r_0 = au
    d.h_0 = 0.01 * au
    d.p = -1.
    d.beta = 1.25
    d.dust = dust

    m.set_spherical_polar_grid_auto(2, 2, 2)

    m.set_n_photons(initial=100, imaging=100)

    m.write(random_filename())


def test_nogrid_invalid():

    m = AnalyticalYSOModel()
    with pytest.raises(Exception) as e:
        m.write(random_filename())
    assert e.value.message == 'The coordinate grid needs to be defined before calling AnalyticalModelYSO.write(...)'


def test_nostar_invalid():

    m = AnalyticalYSOModel()
    with pytest.raises(Exception) as e:
        m.set_spherical_polar_grid_auto(1, 1, 1)
    assert e.value.message == 'The central source radius need to be defined before the grid can be set up'
