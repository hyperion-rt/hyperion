from __future__ import print_function, division

import shutil

import numpy as np
import pytest

import h5py

from .. import Model
from ...grid import AMRGrid
from ...util.functions import random_id
from .test_helpers import get_test_dust, get_realistic_test_dust


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


# Default bin edges for the tests, spanning the full frequency range of both
# test dusts so that the spectrum captures all of the absorbed energy.
_DEFAULT_EDGES = np.logspace(6., 18., 13)


def _read_output_bin_edges(filename):
    """Return the bin edges dataset written to the output (the last matching
    dataset, which is the output one rather than the copied input)."""
    return _read_dataset(filename, '/specific_energy_spectrum_bin_edges')


def _assert_spectrum_sums_to_specific_energy(filename, exclude_empty_spectrum=False):
    # Summed over frequency, specific_energy_spectrum must reproduce specific_energy.
    # Cells with no absorption are clamped up to the minimum specific energy
    # (which only affects specific_energy, not the unclamped specific_energy_spectrum),
    # so we compare only cells that were genuinely heated above that floor.
    # With exclude_empty_spectrum, cells whose spectrum is identically zero are
    # also skipped: the PDA can heat cells that received no photons at all, and
    # for those there is no Monte-Carlo spectral shape to rescale.
    se = _read_dataset(filename, '/specific_energy')
    se_nu = _read_dataset(filename, '/specific_energy_spectrum')
    assert se_nu is not None
    nu_sum = se_nu.sum(axis=0)
    floor = se.min()
    heated = se > floor * (1. + 1.e-6)
    if exclude_empty_spectrum:
        heated &= nu_sum > 0.
    assert np.count_nonzero(heated) >= 1
    np.testing.assert_allclose(nu_sum[heated], se[heated], rtol=1.e-6)


def _cartesian_model(output_specific_energy_spectrum='last', n_cells=2, n_photons=100000,
                     density=1.e-18, bin_edges=None):
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
    m.conf.output.output_specific_energy_spectrum = output_specific_energy_spectrum
    if output_specific_energy_spectrum != 'none':
        m.set_specific_energy_spectrum_bins(_DEFAULT_EDGES if bin_edges is None else bin_edges)
    return m


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_does_not_change_specific_energy(tmpdir):
    # Computing the frequency-resolved specific energy must not perturb the
    # ordinary specific-energy (temperature) calculation. This guards against the
    # regression where specific_energy_sum accumulation was gated on the option.
    se = {}
    for output in ('none', 'last'):
        m = _cartesian_model(output_specific_energy_spectrum=output)
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        se[output] = _read_dataset(out.filename, '/specific_energy')
    # Same seed, and the extra code path does not touch the random stream, so
    # the specific energy should be identical with and without the option.
    np.testing.assert_allclose(se['last'], se['none'], rtol=1.e-10)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_sums_to_specific_energy(tmpdir):
    # specific_energy_spectrum, summed over frequency, must reproduce specific_energy.
    m = _cartesian_model('last')
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_sums_to_specific_energy_with_pda(tmpdir):
    # The partial diffusion approximation overwrites specific_energy in cells
    # with few photon deposits, so the frequency-resolved values in those cells
    # have to be rescaled to match, rather than keeping the noisy Monte-Carlo
    # values. Few photons over many cells guarantees PDA cells exist.
    m = _cartesian_model('last', n_cells=8, n_photons=2000)
    m.set_pda(True)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename, exclude_empty_spectrum=True)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_bin_edges_written(tmpdir):
    # The bin edges are a per-frequency (1-D) quantity, not per-cell (this
    # guards against the regression where the frequency metadata was written
    # as a grid array, which segfaulted whenever n_nu was smaller than
    # n_cells). The edges written to the output must match the ones given.
    m = _cartesian_model('last', n_cells=8)  # 512 cells
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    edges = _read_output_bin_edges(out.filename)
    assert edges is not None
    np.testing.assert_allclose(edges, _DEFAULT_EDGES, rtol=1.e-6)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_dustless_model_runs(tmpdir):
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
def test_specific_energy_spectrum_amr(tmpdir):
    # specific_energy_spectrum must work on AMR grids: it should run and, summed over
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
    m.conf.output.output_specific_energy_spectrum = 'last'
    m.set_specific_energy_spectrum_bins(_DEFAULT_EDGES)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_custom_bin_edges(tmpdir):
    # A user-specified frequency grid given as explicit bin edges: n + 1 edges
    # should produce n bins, the edges should be written to the output, and
    # energy must be conserved when summing over frequency as long as the edges
    # span the frequency range over which energy is absorbed.
    bin_edges = np.logspace(6., 18., 13)
    m = _cartesian_model('last', density=1.e-16, bin_edges=bin_edges)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    edges_out = _read_output_bin_edges(out.filename)
    np.testing.assert_allclose(edges_out, bin_edges, rtol=1.e-6)
    se_nu = _read_dataset(out.filename, '/specific_energy_spectrum')
    assert se_nu.shape[0] == len(bin_edges) - 1
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_windowed_bin_edges(tmpdir):
    # With explicit bin edges, energy absorbed from photons with frequencies
    # outside the outer edges is not counted, so with edges covering only part
    # of the frequency range over which energy is absorbed, the sum over the
    # bins should recover only part of the scalar specific energy.
    m_full = _cartesian_model('last', density=1.e-16, bin_edges=np.logspace(5., 20., 16))
    m_window = _cartesian_model('last', density=1.e-16, bin_edges=np.logspace(14.8, 15.2, 5))
    se_sum = {}
    for key, m in (('full', m_full), ('window', m_window)):
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        se = _read_dataset(out.filename, '/specific_energy')
        se_nu = _read_dataset(out.filename, '/specific_energy_spectrum')
        se_sum[key] = se_nu.sum(axis=0)
        # The spectrum can never contain more energy than the scalar total
        assert np.all(se_sum[key] <= se * (1. + 1.e-6))
    # The window contains part, but only part, of the absorbed energy (for a
    # 6000 K source, roughly a quarter of the energy falls in this window)
    assert 0. < se_sum['window'].sum() < 0.9 * se_sum['full'].sum()


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_with_mrw(tmpdir):
    # Energy deposited during modified-random-walk steps must appear in the
    # frequency-resolved spectrum as well as in the scalar specific energy.
    # This guards against the bug where grid_do_mrw updated only
    # specific_energy_sum, so the spectrum was biased low in the optically
    # thick cells where the MRW handles most of the absorption.
    m = Model()
    m.set_cartesian_grid([-1., 0., 1.], [-1., 0., 1.], [-1., 0., 1.])
    # The density is chosen so that the cells are optically thick enough to
    # the local Planck-mean opacity for MRW steps to actually occur.
    m.add_density_grid(np.ones((2, 2, 2)) * 1.e5, get_realistic_test_dust())
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=1000, imaging=0)
    m.set_seed(-12345)
    m.set_mrw(True, gamma=2.)
    m.set_max_interactions(1000000000)
    m.conf.output.output_specific_energy = 'last'
    m.conf.output.output_specific_energy_spectrum = 'last'
    m.set_specific_energy_spectrum_bins(_DEFAULT_EDGES)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_with_sublimation_cap(tmpdir):
    # With a dust type using the 'cap' sublimation mode, cells below the
    # sublimation threshold must keep their frequency-resolved spectrum. This
    # guards against the bug where sublimate_dust reset the spectrum in every
    # cell rather than only in cells above the threshold, wiping the output.
    m = Model()
    m.set_cartesian_grid([-1., 0., 1.], [-1., 0., 1.], [-1., 0., 1.])
    dust = get_test_dust()
    dust.set_sublimation_temperature('cap', temperature=1600.)
    m.add_density_grid(np.ones((2, 2, 2)) * 1.e-16, dust)
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=100000, imaging=0)
    m.set_seed(-12345)
    m.conf.output.output_specific_energy = 'last'
    m.conf.output.output_specific_energy_spectrum = 'last'
    m.set_specific_energy_spectrum_bins(_DEFAULT_EDGES)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


def test_specific_energy_spectrum_bins_roundtrip(tmpdir):
    # The frequency bin edges should survive a write/read round-trip.
    from ...conf.conf_files import RunConf
    edges = np.logspace(10., 16., 8)
    conf = RunConf()
    conf.set_specific_energy_spectrum_bins(edges)
    with h5py.File(tmpdir.join('conf.h5').strpath, 'w') as f:
        conf._write_specific_energy_spectrum_bins(f)
        conf2 = RunConf()
        conf2._read_specific_energy_spectrum_bins(f)
    np.testing.assert_allclose(conf2.specific_energy_spectrum_bin_edges, edges)


def test_specific_energy_spectrum_bins_validation():
    from ...conf.conf_files import RunConf
    conf = RunConf()
    with pytest.raises(TypeError):
        conf.set_specific_energy_spectrum_bins()
    with pytest.raises(ValueError, match='at least two'):
        conf.set_specific_energy_spectrum_bins([1.e10])
    with pytest.raises(ValueError, match='positive'):
        conf.set_specific_energy_spectrum_bins([-1.e10, 1.e12])
    with pytest.raises(ValueError, match='strictly increasing'):
        conf.set_specific_energy_spectrum_bins([1.e10, 1.e16, 1.e12])
    with pytest.raises(ValueError, match='1-d'):
        conf.set_specific_energy_spectrum_bins([[1.e10, 1.e12], [1.e14, 1.e16]])
    # Default: no bins set
    assert conf.specific_energy_spectrum_bin_edges is None


def test_specific_energy_spectrum_requires_bins(tmpdir):
    # Enabling the frequency-resolved specific energy output without setting
    # the bins must raise a clear error at write time.
    m = _cartesian_model('none')
    m.conf.output.output_specific_energy_spectrum = 'last'
    with pytest.raises(ValueError, match='frequency bins have not been set'):
        m.write(tmpdir.join(random_id()).strpath)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_voronoi(tmpdir):
    # specific_energy_spectrum must also work on (unstructured) Voronoi grids, where
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
    m.conf.output.output_specific_energy_spectrum = 'last'
    m.set_specific_energy_spectrum_bins(_DEFAULT_EDGES)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_get_quantities(tmpdir):
    # specific_energy_spectrum should be retrievable through get_quantities, the
    # per-frequency frequencies metadata should not pollute the grid quantities,
    # and the retrieved array should still conserve energy.
    m = _cartesian_model('last', density=1.e-16)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    g = out.get_quantities()
    assert 'specific_energy_spectrum' in g.quantities
    assert 'specific_energy_spectrum_bin_edges' not in g.quantities
    se = np.array(g.quantities['specific_energy'])       # (dust, nz, ny, nx)
    se_nu = np.array(g.quantities['specific_energy_spectrum'])  # (nu, dust, nz, ny, nx)
    assert se_nu.shape[0] == len(_DEFAULT_EDGES) - 1  # number of bins
    nu_sum = se_nu.sum(axis=0)
    heated = se > se.min() * (1. + 1.e-6)
    np.testing.assert_allclose(nu_sum[heated], se[heated], rtol=1.e-6)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_yt_export_skips_frequency_resolved(tmpdir):
    # The frequency-resolved specific_energy_spectrum is not a per-cell scalar, so it
    # must be skipped (not exported as a malformed field) when converting to yt,
    # while the ordinary scalar quantities are still exported.
    pytest.importorskip("yt")
    m = _cartesian_model('last', density=1.e-16)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    g = out.get_quantities()
    assert 'specific_energy_spectrum' in g.quantities
    ds = g.to_yt()
    fields = [str(f) for f in ds.field_list]
    assert not any('specific_energy_spectrum' in f for f in fields)
    assert any(f.endswith("'density')") for f in fields)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_mpi_matches_serial(tmpdir):
    # The saved spectrum must not depend on the number of MPI processes. We
    # compare the total (which is conserved, hence robust to Monte Carlo noise)
    # between a serial run and a 2-process MPI run. This is skipped unless both an
    # MPI launcher and the MPI build of the binary are available (e.g. CI only
    # builds the serial binaries).
    if shutil.which('mpirun') is None and shutil.which('mpiexec') is None:
        pytest.skip("no MPI launcher available")
    if shutil.which('hyperion_car_mpi') is None:
        pytest.skip("MPI build of the Hyperion binaries not available")

    m = _cartesian_model('last', n_photons=200000)
    m.write(tmpdir.join(random_id()).strpath)

    out_serial = m.run(tmpdir.join(random_id()).strpath)
    out_mpi = m.run(tmpdir.join(random_id()).strpath, mpi=True, n_processes=2)

    total_serial = np.nansum(_read_dataset(out_serial.filename, '/specific_energy_spectrum'))
    total_mpi = np.nansum(_read_dataset(out_mpi.filename, '/specific_energy_spectrum'))
    np.testing.assert_allclose(total_mpi, total_serial, rtol=2.e-2)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_sums_to_specific_energy_with_mrw(tmpdir):
    # In cells heated via the Modified Random Walk, the deposited energy is
    # distributed over the frequency bins according to the local emissivity
    # (the radiation field is Planckian in the diffusion regime), so the sum
    # over frequency must still reproduce the scalar specific energy. The
    # high density here makes the model optically thick enough for the MRW
    # to handle most of the energy deposition.
    m = _cartesian_model('last', n_photons=2000, density=8.)
    m.set_mrw(True, gamma=2.)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_spectrum_sums_to_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_specific_energy_spectrum_is_passive_with_mrw(tmpdir):
    # The MRW spectrum deposit distributes energy deterministically instead
    # of sampling a frequency, so enabling the spectrum must not consume
    # random numbers and therefore must not change the specific energy at
    # all, even in MRW-dominated models.
    se = {}
    for output in ('none', 'last'):
        m = _cartesian_model(output_specific_energy_spectrum=output,
                             n_photons=2000, density=8.)
        m.set_mrw(True, gamma=2.)
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        se[output] = _read_dataset(out.filename, '/specific_energy')
    np.testing.assert_array_equal(se['last'], se['none'])
