# Test that parameters are preserved when written out and read in again

from __future__ import print_function, division

from itertools import product

from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import OutputConf, RunConf, ImageConf, BinnedImageConf, PeeledImageConf
from ...util.functions import virtual_file


@pytest.mark.parametrize(('attribute', 'value'),
                         list(product(['output_density', 'output_density_diff',
                         'output_specific_energy', 'output_n_photons'],
                         ['none', 'last', 'all'])))
def test_io_output_conf(attribute, value):
    o1 = OutputConf()
    setattr(o1, attribute, value)
    v = virtual_file()
    o1.write(v)
    o2 = OutputConf.read(v)
    assert getattr(o2, attribute) == value


def test_io_image_conf():
    i1 = ImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_output_bytes(4)
    i1.set_track_origin('basic')
    i1.set_uncertainties(True)
    v = virtual_file()
    i1.write(v)
    i2 = ImageConf.read(v)
    assert i2.n_x == i1.n_x
    assert i2.n_y == i1.n_y
    assert i2.xmin == i1.xmin
    assert i2.xmax == i1.xmax
    assert i2.ymin == i1.ymin
    assert i2.ymax == i1.ymax
    assert i2.n_ap == i1.n_ap
    assert i2.ap_min == i1.ap_min
    assert i2.ap_max == i1.ap_max
    assert i2.n_wav == i1.n_wav
    assert i2.wav_min == i1.wav_min
    assert i2.wav_max == i1.wav_max
    assert i2.io_bytes == i1.io_bytes
    assert i2.track_origin == i1.track_origin
    assert i2.uncertainties == i1.uncertainties


def test_io_binned_image_conf():
    i1 = BinnedImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_bins(76, 22)
    v = virtual_file()
    i1.write(v)
    i2 = BinnedImageConf.read(v)
    assert i2.n_theta == i1.n_theta
    assert i2.n_phi == i1.n_phi


def test_io_peeled_image_conf():
    i1 = PeeledImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, 4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_angles([1., 2., 3], [4., 5., 6.])
    i1.set_peeloff_origin([2.2, 3.3, 7.6])
    i1.set_ignore_optical_depth(True)
    i1.set_depth(-1.7, 6.2)
    v = virtual_file()
    i1.write(v)
    i2 = PeeledImageConf.read(v)
    for i in range(len(i2.viewing_angles)):
        assert i2.viewing_angles[i][0] == i1.viewing_angles[i][0]
        assert i2.viewing_angles[i][1] == i1.viewing_angles[i][1]
    assert_equal(i2.peeloff_origin, i1.peeloff_origin)
    assert i2.ignore_optical_depth == i1.ignore_optical_depth
    assert i2.d_min == i1.d_min
    assert i2.d_max == i1.d_max


def test_io_peeled_image_conf_inside():
    i1 = PeeledImageConf()
    i1.set_image_size(33, 42)
    i1.set_image_limits(3.2, -4.4, 5.2, 9.9)
    i1.set_aperture_range(6, 1.2, 8.8)
    i1.set_wavelength_range(9, 2.2, 7.4)
    i1.set_viewing_angles([1., 2., 3], [4., 5., 6.])
    i1.set_inside_observer([7., 8., 9.])
    i1.set_ignore_optical_depth(True)
    i1.set_depth(1.7, 6.2)
    v = virtual_file()
    i1.write(v)
    i2 = PeeledImageConf.read(v)
    for i in range(len(i2.viewing_angles)):
        assert i2.viewing_angles[i][0] == i1.viewing_angles[i][0]
        assert i2.viewing_angles[i][1] == i1.viewing_angles[i][1]
    assert_equal(i2.inside_observer, i1.inside_observer)
    assert i2.ignore_optical_depth == i1.ignore_optical_depth
    assert i2.d_min == i1.d_min
    assert i2.d_max == i1.d_max


# RUNTIME CONFIGURATION

@pytest.mark.parametrize(('value'), [0.001, 0.1, 1.])
def test_io_run_conf_propagation_check_frequency(value):
    r1 = RunConf()
    r1.set_propagation_check_frequency(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2._frequency == r1._frequency


@pytest.mark.parametrize(('value'), [-1234, -6663121])
def test_io_run_conf_seed(value):
    r1 = RunConf()
    r1.set_seed(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2._seed == r1._seed


@pytest.mark.parametrize(('value'), [2, 5, 102])
def test_io_run_conf_n_initial_iterations(value):
    r1 = RunConf()
    r1.set_n_initial_iterations(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.n_iter == r1.n_iter


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_raytracing(value):
    r1 = RunConf()
    r1.set_raytracing(value)
    if value:
        r1.set_n_photons(1, 2, raytracing_sources=3, raytracing_dust=4)
    else:
        r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.raytracing == r1.raytracing


def test_io_run_conf_n_photons_plain():
    r1 = RunConf()
    r1.set_n_photons(initial=1, imaging=2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    for key in r1.n_photons:
        assert r2.n_photons[key] == r1.n_photons[key]


def test_io_run_conf_n_photons_raytracing():
    r1 = RunConf()
    r1.set_raytracing(True)
    r1.set_n_photons(initial=1, imaging=2, raytracing_sources=3,
                     raytracing_dust=4)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    for key in r1.n_photons:
        assert r2.n_photons[key] == r1.n_photons[key]


def test_io_run_conf_n_photons_monochromatic():
    r1 = RunConf()
    r1._monochromatic = True
    r1.set_n_photons(initial=1, imaging_sources=3, imaging_dust=4)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2._monochromatic = True
    r2.read_run_conf(v)
    for key in r1.n_photons:
        assert r2.n_photons[key] == r1.n_photons[key]


@pytest.mark.parametrize(('value'), [33, 5283])
def test_io_run_conf_max_interactions(value):
    r1 = RunConf()
    r1.set_max_interactions(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.n_inter_max == r1.n_inter_max


@pytest.mark.parametrize(('value'), [77, 1244])
def test_io_run_conf_max_reabsorptions(value):
    r1 = RunConf()
    r1.set_max_reabsorptions(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.n_reabs_max == r1.n_reabs_max


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_pda(value):
    r1 = RunConf()
    r1.set_pda(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.pda == r1.pda


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_mrw(value):
    r1 = RunConf()
    r1.set_mrw(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.mrw == r1.mrw


@pytest.mark.parametrize(('value'), [False, True])
def test_io_run_conf_convergence(value):
    r1 = RunConf()
    if value:
        r1.set_convergence(value, percentile=12., absolute=34., relative=56.)
    else:
        r1.set_convergence(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.check_convergence == r1.check_convergence
    assert r2.convergence_percentile == r1.convergence_percentile
    assert r2.convergence_absolute == r1.convergence_absolute
    assert r2.convergence_relative == r1.convergence_relative


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_kill_on_absorb(value):
    r1 = RunConf()
    r1.set_kill_on_absorb(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.kill_on_absorb == r1.kill_on_absorb


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_kill_on_scatter(value):
    r1 = RunConf()
    r1.set_kill_on_scatter(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.kill_on_scatter == r1.kill_on_scatter


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_forced_first_scattering(value):
    r1 = RunConf()
    r1.set_forced_first_scattering(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.forced_first_scattering == r1.forced_first_scattering


@pytest.mark.parametrize(('value'), [4, 8])
def test_io_run_conf_output_bytes(value):
    r1 = RunConf()
    r1.set_output_bytes(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.physics_io_bytes == r1.physics_io_bytes


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_sample_sources_evenly(value):
    r1 = RunConf()
    r1.set_sample_sources_evenly(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.sample_sources_evenly == r1.sample_sources_evenly


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_enforce_energy_range(value):
    r1 = RunConf()
    r1.set_enforce_energy_range(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.enforce_energy_range == r1.enforce_energy_range


@pytest.mark.parametrize(('value'), [True, False])
def test_io_run_conf_copy_input(value):
    r1 = RunConf()
    r1.set_copy_input(value)
    r1.set_n_photons(1, 2)
    v = virtual_file()
    r1.write_run_conf(v)
    r2 = RunConf()
    r2.read_run_conf(v)
    assert r2.copy_input == r1.copy_input
