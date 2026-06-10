from __future__ import print_function, division

import shutil

import numpy as np
import pytest

import h5py

from .. import Model
from ...grid import AMRGrid
from ...util.functions import random_id
from .test_helpers import get_test_dust


def _read_dataset(filename, suffix):
    """Return the last grid dataset whose HDF5 path ends with ``suffix``."""
    with h5py.File(filename, 'r') as f:
        matches = []
        f.visititems(lambda name, obj: matches.append(name)
                     if name.endswith(suffix) and isinstance(obj, h5py.Dataset)
                     else None)
        if not matches:
            return None
        return np.asarray(f[matches[-1]][()], dtype=float)


def _assert_nu_sums_to_specific_energy(filename):
    # Summed over frequency, specific_energy_nu must reproduce specific_energy.
    # Cells with no absorption are clamped up to the minimum specific energy
    # (which only affects specific_energy, not the unclamped specific_energy_nu),
    # so we compare only cells that were genuinely heated above that floor.
    se = _read_dataset(filename, '/specific_energy')
    se_nu = _read_dataset(filename, '/specific_energy_nu')
    assert se_nu is not None
    nu_sum = se_nu.sum(axis=0)
    floor = se.min()
    heated = se > floor * (1. + 1.e-6)
    assert np.count_nonzero(heated) >= 1
    np.testing.assert_allclose(nu_sum[heated], se[heated], rtol=1.e-6)


def _cartesian_model(output_specific_energy_nu='last', n_cells=2, n_photons=100000,
                     density=1.e-18, frequencies=None):
    m = Model()
    edges = np.linspace(-1., 1., n_cells + 1)
    m.set_cartesian_grid(edges, edges, edges)
    dust = get_test_dust()
    m.add_density_grid(np.ones((n_cells, n_cells, n_cells)) * density, dust)
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=n_photons, imaging=0)
    m.set_seed(-12345)
    m.conf.output.output_specific_energy = 'last'
    m.conf.output.output_specific_energy_nu = output_specific_energy_nu
    if frequencies is not None:
        m.set_specific_energy_nu_frequencies(frequencies)
    return m


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_does_not_change_specific_energy(tmpdir):
    # Computing the frequency-resolved specific energy must not perturb the
    # ordinary specific-energy (temperature) calculation. This guards against the
    # regression where specific_energy_sum accumulation was gated on the option.
    se = {}
    for output in ('none', 'last'):
        m = _cartesian_model(output_specific_energy_nu=output)
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        se[output] = _read_dataset(out.filename, '/specific_energy')
    # Same seed, and the extra code path does not touch the random stream, so
    # the specific energy should be identical with and without the option.
    np.testing.assert_allclose(se['last'], se['none'], rtol=1.e-10)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_sums_to_specific_energy(tmpdir):
    # specific_energy_nu, summed over frequency, must reproduce specific_energy.
    m = _cartesian_model('last')
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_nu_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_frequencies_written(tmpdir):
    # specific_energy_nu_frequencies is a per-frequency (1-D) quantity, not
    # per-cell. This guards against the regression where it was written as a grid
    # array, which segfaulted whenever n_nu was smaller than n_cells.
    m = _cartesian_model('last', n_cells=8)  # 512 cells, dust has 2 frequencies
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    bins = _read_dataset(out.filename, '/specific_energy_nu_frequencies')
    assert bins is not None
    assert bins.shape == (2,)  # == number of dust frequencies, not n_cells


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_dustless_model_runs(tmpdir):
    # A model with no dust must still run (and image) without crashing. This
    # guards against the regression where the instrumentation dereferenced the
    # first dust type unconditionally, segfaulting dustless models.
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([45.], [45.])
    i.set_wavelength_range(5, 0.1, 100.)
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=10000)
    m.write(tmpdir.join(random_id()).strpath)
    # Should complete without raising (previously segfaulted in the final iteration).
    m.run(tmpdir.join(random_id()).strpath)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_amr(tmpdir):
    # specific_energy_nu must work on AMR grids: it should run and, summed over
    # frequency, reproduce specific_energy in every cell.
    m = Model()
    amr = AMRGrid()
    level = amr.add_level()
    grid = level.add_grid()
    grid.xmin, grid.xmax = -1., 1.
    grid.ymin, grid.ymax = -1., 1.
    grid.zmin, grid.zmax = -1., 1.
    grid.nx, grid.ny, grid.nz = 4, 4, 4
    grid.quantities['density'] = [np.ones((4, 4, 4)) * 1.e-16]
    m.set_amr_grid(amr)
    m.add_density_grid(amr['density'][0], get_test_dust())
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=100000, imaging=0)
    m.set_seed(-9)
    m.conf.output.output_specific_energy = 'last'
    m.conf.output.output_specific_energy_nu = 'last'
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_nu_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_custom_frequency_grid(tmpdir):
    # A user-specified frequency grid (independent of the dust frequency grid)
    # should be used for the output bins, and energy must still be conserved
    # when summing over frequency.
    frequencies = np.logspace(10., 16., 12)
    m = _cartesian_model('last', density=1.e-16, frequencies=frequencies)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    bins = _read_dataset(out.filename, '/specific_energy_nu_frequencies')
    np.testing.assert_allclose(bins, frequencies, rtol=1.e-6)
    se_nu = _read_dataset(out.filename, '/specific_energy_nu')
    assert se_nu.shape[0] == len(frequencies)
    _assert_nu_sums_to_specific_energy(out.filename)


def test_specific_energy_nu_frequencies_roundtrip(tmpdir):
    # The custom frequency grid should survive a write/read round-trip.
    from ...conf.conf_files import RunConf
    frequencies = np.logspace(10., 16., 8)
    conf = RunConf()
    conf.set_specific_energy_nu_frequencies(frequencies)
    with h5py.File(tmpdir.join('conf.h5').strpath, 'w') as f:
        conf._write_specific_energy_nu_frequencies(f)
        conf2 = RunConf()
        conf2._read_specific_energy_nu_frequencies(f)
    np.testing.assert_allclose(conf2.specific_energy_nu_frequencies, frequencies)


def test_specific_energy_nu_frequencies_validation():
    from ...conf.conf_files import RunConf
    conf = RunConf()
    with pytest.raises(ValueError):
        conf.set_specific_energy_nu_frequencies([[1.e10, 1.e12], [1.e14, 1.e16]])
    with pytest.raises(ValueError):
        conf.set_specific_energy_nu_frequencies([1.e10, -1.e12])
    # Default: no custom grid set
    assert conf.specific_energy_nu_frequencies is None


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_voronoi(tmpdir):
    # specific_energy_nu must also work on (unstructured) Voronoi grids, where
    # cells are a flat list rather than a structured grid.
    rng = np.random.RandomState(12345)
    n = 30
    x = rng.uniform(-1., 1., n)
    y = rng.uniform(-1., 1., n)
    z = rng.uniform(-1., 1., n)
    m = Model()
    m.set_voronoi_grid(x, y, z)
    m.add_density_grid(np.ones(n) * 1.e-16, get_test_dust())
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=100000, imaging=0)
    m.set_seed(-12345)
    m.conf.output.output_specific_energy = 'last'
    m.conf.output.output_specific_energy_nu = 'last'
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_nu_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_get_quantities(tmpdir):
    # specific_energy_nu should be retrievable through get_quantities, the
    # per-frequency frequencies metadata should not pollute the grid quantities,
    # and the retrieved array should still conserve energy.
    m = _cartesian_model('last', density=1.e-16)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    g = out.get_quantities()
    assert 'specific_energy_nu' in g.quantities
    assert 'specific_energy_nu_frequencies' not in g.quantities
    se = np.array(g.quantities['specific_energy'])       # (dust, nz, ny, nx)
    se_nu = np.array(g.quantities['specific_energy_nu'])  # (nu, dust, nz, ny, nx)
    assert se_nu.shape[0] == 2  # number of dust frequencies
    nu_sum = se_nu.sum(axis=0)
    heated = se > se.min() * (1. + 1.e-6)
    np.testing.assert_allclose(nu_sum[heated], se[heated], rtol=1.e-6)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_yt_export_skips_frequency_resolved(tmpdir):
    # The frequency-resolved specific_energy_nu is not a per-cell scalar, so it
    # must be skipped (not exported as a malformed field) when converting to yt,
    # while the ordinary scalar quantities are still exported.
    pytest.importorskip("yt")
    m = _cartesian_model('last', density=1.e-16)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    g = out.get_quantities()
    assert 'specific_energy_nu' in g.quantities
    ds = g.to_yt()
    fields = [str(f) for f in ds.field_list]
    assert not any('specific_energy_nu' in f for f in fields)
    assert any(f.endswith("'density')") for f in fields)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_nu_mpi_matches_serial(tmpdir):
    # The saved spectrum must not depend on the number of MPI processes. We
    # compare the total (which is conserved, hence robust to Monte Carlo noise)
    # between a serial run and a 2-process MPI run.
    if shutil.which('mpirun') is None and shutil.which('mpiexec') is None:
        pytest.skip("no MPI launcher available")

    m = _cartesian_model('last', n_photons=200000)
    m.write(tmpdir.join(random_id()).strpath)

    out_serial = m.run(tmpdir.join(random_id()).strpath)
    out_mpi = m.run(tmpdir.join(random_id()).strpath, mpi=True, n_processes=2)

    total_serial = np.nansum(_read_dataset(out_serial.filename, '/specific_energy_nu'))
    total_mpi = np.nansum(_read_dataset(out_mpi.filename, '/specific_energy_nu'))
    np.testing.assert_allclose(total_mpi, total_serial, rtol=2.e-2)
