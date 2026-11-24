import shutil

import numpy as np
import pytest

from .. import Model
from .test_helpers import random_id
from ...dust import IsotropicDust


def scattering_model(seed):

    # Point source at the centre of an optically thin sphere of mostly
    # scattering dust, so that the scattered flux depends on the random
    # positions and directions of the photons

    m = Model()
    m.set_spherical_polar_grid(np.linspace(0., 1., 21), np.linspace(0., np.pi, 11), [0., 2. * np.pi])
    dust = IsotropicDust([3.e9, 3.e16], [0.9, 0.9], [1., 1.])
    dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)
    m.add_density_grid(np.ones((1, 10, 20)) * 0.3, dust)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([10., 50., 90.], [0., 0., 0.])
    i.set_track_origin('basic')
    i.set_uncertainties(True)

    m.set_n_initial_iterations(0)
    m.set_monochromatic(True, wavelengths=[1.])
    m.set_n_photons(imaging_sources=40000, imaging_dust=0)
    m.set_seed(seed)
    return m


@pytest.mark.requires_hyperion_binaries
def test_mpi_small_seed(tmpdir):

    # Each MPI process needs a different random number sequence. Previously,
    # process k used seed + k, so that with a seed of -2 and 4 processes,
    # one process had a seed of 0, which gives a degenerate sequence of random
    # numbers, and another had a seed of 1, which gives the same sequence as
    # -1 (since the absolute value of the seed is used). The scattered flux
    # was then wrong. This should now agree with a serial run.

    if shutil.which('mpirun') is None and shutil.which('mpiexec') is None:
        pytest.skip("no MPI launcher available")
    if shutil.which('hyperion_sph_mpi') is None:
        pytest.skip("MPI build of the Hyperion binaries not available")

    m = scattering_model(-2)
    m.write(tmpdir.join(random_id()).strpath)
    serial = m.run(tmpdir.join(random_id()).strpath).get_sed(component='source_scat', uncertainties=True)
    mpi = m.run(tmpdir.join(random_id()).strpath, mpi=True, n_processes=4).get_sed(component='source_scat', uncertainties=True)

    diff = np.abs(mpi.val - serial.val) / np.hypot(mpi.unc, serial.unc)
    assert np.all(diff < 5.)


@pytest.mark.parametrize('seed', [-65543, -131086])
def test_degenerate_seed(seed):
    m = Model()
    with pytest.raises(ValueError) as exc:
        m.set_seed(seed)
    assert exc.value.args[0] == "seed should not be a multiple of 65543, which gives a degenerate sequence of random numbers"


@pytest.mark.parametrize('seed', [-1, -2, -65542, -65544, -124902])
def test_valid_seed(seed):
    m = Model()
    m.set_seed(seed)
