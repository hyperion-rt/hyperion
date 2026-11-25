import os

import pytest
import numpy as np

from ..functions import B_nu, dB_nu_dT
from ..integrate import integrate_loglog
from ..constants import sigma, k, c


def test_b_nu():

    nu = np.logspace(-20, 20., 10000)

    for T in [10, 100, 1000, 10000]:

        # Compute planck function
        b = B_nu(nu, T)

        # Check that the intergral is correct
        total = integrate_loglog(nu, b)
        np.testing.assert_allclose(total, sigma * T ** 4 / np.pi, rtol=1e-4)

        # Check that we reach the rayleigh-jeans limit at low frequencies
        rj = 2. * nu ** 2 * k * T / c**2
        np.testing.assert_allclose(b[nu < 1e-10], rj[nu < 1e-10], rtol=1.e-8)


def test_db_nu_dt():

    nu = np.logspace(-20, 20., 10000)

    for T in [10, 100, 1000, 10000]:

        # Compute exact planck function derivative
        db = dB_nu_dT(nu, T)

        # Compute numerical planck function derivative
        dT = T / 1e6
        b1 = B_nu(nu, T - dT)
        b2 = B_nu(nu, T + dT)
        db_num = 0.5 * (b2 - b1) / dT

        # Check that the two are the same
        np.testing.assert_allclose(db, db_num, rtol=1.e-2)


def test_freezable_attributes_per_instance():
    # Regression: FreezableClass used to store _attributes on the class, so
    # attributes registered on one instance leaked into every other instance
    # (and grew an unbounded shared list). They must be per-instance.
    from ..functions import FreezableClass

    class Foo(FreezableClass):
        def __init__(self):
            FreezableClass.__init__(self)

    a = Foo()
    a.x = 1        # registers 'x' on a only
    a._freeze()

    b = Foo()
    b._freeze()
    # 'x' belongs to a, not b; with the shared-list bug b would allow it
    with pytest.raises(AttributeError):
        b.x = 2

    # a can still change its own registered attribute
    a.x = 5
    assert a.x == 5


@pytest.mark.parametrize('absolute_paths', [False, True])
def test_link_or_copy_link(tmpdir, absolute_paths):

    import h5py
    from ..functions import link_or_copy

    target = tmpdir.mkdir('a').join('target.hdf5').strpath
    with h5py.File(target, 'w') as f:
        f['x'] = np.arange(5.)

    filename = tmpdir.mkdir('b').join('new.hdf5').strpath
    with h5py.File(filename, 'w') as f:
        link_or_copy(f, 'x', h5py.ExternalLink(target, '/x'), False, absolute_paths=absolute_paths)

    with h5py.File(filename, 'r') as f:
        link = f.get('x', getlink=True)
        assert isinstance(link, h5py.ExternalLink)
        assert link.filename == (target if absolute_paths else os.path.join('..', 'a', 'target.hdf5'))
        np.testing.assert_equal(f['x'][()], np.arange(5.))


def test_link_or_copy_link_not_opened(tmpdir):

    # The target of the link should not be opened when creating the link, since
    # it would be opened in read/write mode (and therefore locked) when
    # following the link from a file open for writing, which fails if another
    # process has the target open, and can modify the target.

    import h5py
    from ..functions import link_or_copy

    filename = tmpdir.join('new.hdf5').strpath
    with h5py.File(filename, 'w') as f:
        link_or_copy(f, 'x', h5py.ExternalLink(tmpdir.join('missing.hdf5').strpath, '/x'), False)
        assert isinstance(f.get('x', getlink=True), h5py.ExternalLink)
