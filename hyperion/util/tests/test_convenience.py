from __future__ import print_function, division

from ..convenience import OptThinRadius


def test_optthin_radius_min_stored():
    r = OptThinRadius(1600., min=3.)
    assert r.temperature == 1600.
    assert r.min == 3.


def test_optthin_radius_min_propagated():
    r = OptThinRadius(1600., min=3.)
    assert (r * 2.).min == 3.
    assert (2. * r).min == 3.
