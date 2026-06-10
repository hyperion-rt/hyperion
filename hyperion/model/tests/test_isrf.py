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


def _assert_isrf_reconstructs_specific_energy(filename):
    # Summed over frequency, the ISRF must reproduce the specific energy. Cells
    # with no absorption are clamped up to the minimum specific energy (which
    # only affects specific_energy, not the unclamped specific_energy_nu), so we
    # compare only cells that were genuinely heated above that floor.
    se = _read_dataset(filename, '/specific_energy')
    se_nu = _read_dataset(filename, '/specific_energy_nu')
    assert se_nu is not None
    nu_sum = se_nu.sum(axis=0)
    floor = se.min()
    heated = se > floor * (1. + 1.e-6)
    assert np.count_nonzero(heated) >= 1
    np.testing.assert_allclose(nu_sum[heated], se[heated], rtol=1.e-6)


def _cartesian_isrf_model(compute_isrf, n_cells=2, n_photons=100000, density=1.e-18):
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
    m.compute_isrf(compute_isrf)
    return m


@pytest.mark.requires_hyperion_binaries
def test_isrf_does_not_change_specific_energy(tmpdir):
    # Computing the ISRF must not perturb the ordinary specific-energy
    # (temperature) calculation. This guards against the regression where the
    # specific_energy_sum accumulation was accidentally gated behind compute_isrf.
    se = {}
    for isrf in (False, True):
        m = _cartesian_isrf_model(isrf)
        m.write(tmpdir.join(random_id()).strpath)
        out = m.run(tmpdir.join(random_id()).strpath)
        se[isrf] = _read_dataset(out.filename, '/specific_energy')
    # Same seed, and the ISRF code path does not touch the random stream, so
    # the specific energy should be identical with and without the ISRF.
    np.testing.assert_allclose(se[True], se[False], rtol=1.e-10)


@pytest.mark.requires_hyperion_binaries
def test_isrf_sums_to_specific_energy(tmpdir):
    # The frequency-resolved ISRF, summed over frequency, must reproduce the
    # total specific energy.
    m = _cartesian_isrf_model(True)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_isrf_reconstructs_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_isrf_frequency_bins_written(tmpdir):
    # ISRF_frequency_bins is a per-frequency (1-D) quantity, not per-cell. This
    # guards against the regression where it was written as a grid array, which
    # segfaulted whenever the number of frequencies was smaller than n_cells.
    m = _cartesian_isrf_model(True, n_cells=8)  # 512 cells, dust has 2 frequencies
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    bins = _read_dataset(out.filename, '/ISRF_frequency_bins')
    assert bins is not None
    assert bins.shape == (2,)  # == number of dust frequencies, not n_cells


@pytest.mark.requires_hyperion_binaries
def test_isrf_dustless_model_runs(tmpdir):
    # A model with no dust must still run (and image) without crashing. This
    # guards against the regression where ISRF instrumentation dereferenced the
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
def test_isrf_amr(tmpdir):
    # The ISRF must work on AMR grids: it should run and, summed over frequency,
    # reproduce the specific energy in every cell.
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
    m.compute_isrf(True)
    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    _assert_isrf_reconstructs_specific_energy(out.filename)


@pytest.mark.requires_hyperion_binaries
def test_isrf_mpi_matches_serial(tmpdir):
    # The saved ISRF must not depend on the number of MPI processes. We compare
    # the total ISRF energy (which is conserved, hence robust to Monte Carlo
    # noise) between a serial run and a 2-process MPI run.
    if shutil.which('mpirun') is None and shutil.which('mpiexec') is None:
        pytest.skip("no MPI launcher available")

    m = _cartesian_isrf_model(True, n_photons=200000)
    m.write(tmpdir.join(random_id()).strpath)

    out_serial = m.run(tmpdir.join(random_id()).strpath)
    out_mpi = m.run(tmpdir.join(random_id()).strpath, mpi=True, n_processes=2)

    total_serial = np.nansum(_read_dataset(out_serial.filename, '/specific_energy_nu'))
    total_mpi = np.nansum(_read_dataset(out_mpi.filename, '/specific_energy_nu'))
    np.testing.assert_allclose(total_mpi, total_serial, rtol=2.e-2)
