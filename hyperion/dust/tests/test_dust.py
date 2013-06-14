import os
import tempfile

from astropy.tests.helper import pytest

from .. import SphericalDust
from ...util.functions import random_id


def test_missing_properties(tmpdir):
    d = SphericalDust()
    with pytest.raises(Exception) as e:
        d.write(tmpdir.join(random_id()).strpath)
    assert e.value.args[0] == "The following attributes of the optical properties have not been set: nu, chi, albedo, mu, P1, P2, P3, P4"
