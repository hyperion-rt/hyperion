import os
import tempfile

import pytest

from .. import SphericalDust
from ...util.functions import random_id


def random_filename():
    return os.path.join(tempfile.mkdtemp(), random_id())


def test_missing_properties():
    d = SphericalDust()
    with pytest.raises(Exception) as e:
        d.write(random_filename())
    assert e.value.args[0] == "The following attributes of the optical properties have not been set: nu, chi, albedo, mu, P1, P2, P3, P4"
