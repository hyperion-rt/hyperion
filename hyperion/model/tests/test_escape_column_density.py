import numpy as np
import pytest

from .. import Model
from .test_helpers import random_id
from ...dust import IsotropicDust

# The models in this file are uniform-density spheres of radius 1 divided into
# radial shells, for which the minimum column density to escape from anywhere
# in a shell is reached by going radially outwards from its outer wall, so
# that it is equal to the density times (1 - r_outer).

N_SHELLS = 50
RHO = 5.
CHI = 1.


def sphere_model(peeloff_tau_max=None, escape=True, emission_max_escape_optical_depth=None,
                 specific_energy=None, dust_photons=False):

    m = Model()
    m.set_spherical_polar_grid(np.linspace(0., 1., N_SHELLS + 1), [0., np.pi], [0., 2. * np.pi])

    dust = IsotropicDust([3.e9, 3.e16], [0., 0.], [CHI, CHI])
    dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)

    density = np.ones((1, 1, N_SHELLS)) * RHO
    if specific_energy is None:
        m.add_density_grid(density, dust)
    else:
        m.add_density_grid(density, dust, specific_energy=density * 0. + specific_energy)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([90.], [0.])
    i.set_track_origin('basic')

    m.set_n_initial_iterations(0)
    m.set_monochromatic(True, wavelengths=[1.],
                        emission_max_escape_optical_depth=emission_max_escape_optical_depth)

    if dust_photons:
        m.set_n_photons(imaging_sources=0, imaging_dust=10000)
    else:
        m.set_n_photons(imaging_sources=1000, imaging_dust=0)

    if peeloff_tau_max is not None:
        m.set_max_peeloff_optical_depth(peeloff_tau_max)

    if escape:
        m.set_compute_escape_column_density(True, rays_per_cell=30)

    return m


def run(tmpdir, m):
    m.write(tmpdir.join(random_id()).strpath)
    return m.run(tmpdir.join(random_id()).strpath)


@pytest.mark.requires_hyperion_binaries
def test_escape_column_density_sphere(tmpdir):

    out = run(tmpdir, sphere_model())

    column = out.get_quantities()['minimum_escape_column_density'][0].array[0, 0, :]

    r_outer = np.linspace(0., 1., N_SHELLS + 1)[1:]
    expected = RHO * (1. - r_outer)

    # The estimate can never be lower than the true minimum, and with 30
    # rays per cell should be very close to it
    assert np.all(column >= expected - 1e-10)
    np.testing.assert_allclose(column, expected, rtol=0.01, atol=1e-10)

    # The outermost shell is on the edge of the grid
    assert column[-1] == 0.



@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize('grid_type', ['car', 'amr'])
def test_escape_column_density_cube(tmpdir, grid_type):

    # A uniform cube of side 2 divided into 8x8x8 cells, for which the minimum
    # column density to escape a cell is the density times the distance from
    # the cell to the nearest face of the cube. This also checks that the
    # minimum escape column density is written and read correctly for grid
    # types other than spherical polar grids.

    from ...grid import AMRGrid

    m = Model()
    if grid_type == 'car':
        w = np.linspace(-1., 1., 9)
        m.set_cartesian_grid(w, w, w)
        m.add_density_grid(np.ones((8, 8, 8)) * RHO, IsotropicDust([3.e9, 3.e16], [0., 0.], [CHI, CHI]))
    else:
        amr = AMRGrid()
        grid = amr.add_level().add_grid()
        grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax = -1., 1., -1., 1., -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 8, 8
        grid.quantities['density'] = np.ones((8, 8, 8)) * RHO
        m.set_amr_grid(amr)
        m.add_density_grid(amr['density'], IsotropicDust([3.e9, 3.e16], [0., 0.], [CHI, CHI]))
    for dust in m.dust:
        dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=0)
    m.set_compute_escape_column_density(True, rays_per_cell=100)

    out = run(tmpdir, m)
    column = out.get_quantities()['minimum_escape_column_density'][0]
    if grid_type == 'car':
        column = column.array
    else:
        column = column.levels[0].grids[0].quantities['minimum_escape_column_density']

    # Cell centres along each axis, and the distance from the centre to the
    # walls (0.125)
    c = np.linspace(-0.875, 0.875, 8)
    z, y, x = np.meshgrid(c, c, c, indexing='ij')
    expected = RHO * (np.minimum.reduce([1. - np.abs(x), 1. - np.abs(y), 1. - np.abs(z)]) - 0.125)

    assert np.all(column >= expected - 1e-10)
    np.testing.assert_allclose(column, expected, rtol=0.02, atol=1e-10)


@pytest.mark.requires_hyperion_binaries
def test_escape_column_density_not_computed(tmpdir):
    out = run(tmpdir, sphere_model(escape=False))
    assert 'minimum_escape_column_density' not in out.get_quantities()


def source_flux(tmpdir, **kwargs):
    out = run(tmpdir, sphere_model(**kwargs))
    return out.get_sed(component='source_emit').val[0, 0, 0]


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize('escape', [False, True])
def test_peeloff_hidden_cells(tmpdir, escape):

    # The source is in the central cell, from which the minimum optical depth
    # to escape is (1 - 1 / N_SHELLS) * RHO * CHI = 4.9, and the optical depth
    # to the observer is RHO * CHI = 5.

    reference = source_flux(tmpdir, escape=False, peeloff_tau_max=1000.) / np.exp(-5.)

    # The optical depth to escape is lower than the maximum, so the photons
    # are peeled off as usual
    flux = source_flux(tmpdir, escape=escape, peeloff_tau_max=10.)
    np.testing.assert_allclose(flux / reference, np.exp(-5.), rtol=1e-10)

    # The optical depth to escape is higher than the maximum. Without the
    # minimum escape column density, the optical depth integration stops once
    # it exceeds the maximum, and with it, the photons are not peeled off at
    # all.
    flux = source_flux(tmpdir, escape=escape, peeloff_tau_max=2.)
    if escape:
        assert flux == 0.
    else:
        assert np.exp(-2.1) <= flux / reference <= np.exp(-2.) * (1 + 1e-10)


def dust_flux(tmpdir, emission_max_escape_optical_depth):
    out = run(tmpdir, sphere_model(dust_photons=True, specific_energy=1.e10,
                                   emission_max_escape_optical_depth=emission_max_escape_optical_depth))
    return out.get_sed(component='dust_emit').val[0, 0, 0]


@pytest.mark.requires_hyperion_binaries
def test_emission_hidden_cells(tmpdir):

    # The maximum optical depth to escape any cell is 4.9, so with a maximum
    # of 10, no cells should be removed, and the results should be identical
    # to the case without any maximum
    flux_all = dust_flux(tmpdir, None)
    assert flux_all > 0.
    assert dust_flux(tmpdir, 10.) == flux_all

    # With a very low maximum, only the outermost shell, from which the
    # optical depth to escape is zero, emits
    flux_outer = dust_flux(tmpdir, 1e-10)
    assert 0. < flux_outer < flux_all


@pytest.mark.parametrize('image_setting', ['inside_observer', 'ignore_optical_depth'])
def test_emission_hidden_cells_invalid_images(tmpdir, image_setting):

    # Cells hidden from observers outside the grid can be visible to inside
    # observers and in images that ignore the optical depth, so skipping their
    # emission is not allowed in these cases

    m = sphere_model(dust_photons=True, emission_max_escape_optical_depth=200.)
    i = m.add_peeled_images(sed=False, image=True)
    i.set_image_size(1, 1)
    if image_setting == 'inside_observer':
        i.set_inside_observer((0., 0., 0.))
        i.set_image_limits(180., -180., -90., 90.)
        message = "with inside observers"
    else:
        i.set_viewing_angles([90.], [0.])
        i.set_image_limits(-1., 1., -1., 1.)
        i.set_ignore_optical_depth(True)
        message = "with peeled images that ignore the optical depth"

    with pytest.raises(ValueError) as exc:
        m.write(tmpdir.join(random_id()).strpath)
    assert exc.value.args[0] == "Cannot use emission_max_escape_optical_depth " + message

    # Without emission_max_escape_optical_depth, the model can be written
    m.set_monochromatic(True, wavelengths=[1.])
    m.write(tmpdir.join(random_id()).strpath)


def test_emission_hidden_cells_requires_escape_column_density(tmpdir):
    m = sphere_model(escape=False, dust_photons=True, emission_max_escape_optical_depth=200.)
    with pytest.raises(ValueError) as exc:
        m.write(tmpdir.join(random_id()).strpath)
    assert exc.value.args[0] == ("emission_max_escape_optical_depth requires the minimum escape column "
                                 "density to be computed - use set_compute_escape_column_density to enable it")
    m.set_compute_escape_column_density(True, rays_per_cell=30)
    m.write(tmpdir.join(random_id()).strpath)


def test_emission_hidden_cells_requires_monochromatic():
    m = Model()
    with pytest.raises(ValueError) as exc:
        m.set_monochromatic(False, emission_max_escape_optical_depth=200.)
    assert exc.value.args[0] == "emission_max_escape_optical_depth can only be set in monochromatic mode"


@pytest.mark.parametrize(('value', 'exception', 'message'), [
    ('a', TypeError, "emission_max_escape_optical_depth should be a scalar value"),
    ([1., 2.], TypeError, "emission_max_escape_optical_depth should be a scalar value"),
    (0., ValueError, "emission_max_escape_optical_depth should be strictly positive"),
    (float('nan'), ValueError, "emission_max_escape_optical_depth should be strictly positive"),
])
def test_emission_hidden_cells_invalid_value(value, exception, message):
    m = Model()
    with pytest.raises(exception) as exc:
        m.set_monochromatic(True, wavelengths=[1.], emission_max_escape_optical_depth=value)
    assert exc.value.args[0] == message


@pytest.mark.parametrize('value', [None, 200.])
def test_emission_hidden_cells_io(value):
    from ...util.functions import virtual_file
    m1 = Model()
    m1.set_monochromatic(True, wavelengths=[1.], emission_max_escape_optical_depth=value)
    v = virtual_file()
    m1._write_monochromatic(v)
    m2 = Model()
    m2.set_monochromatic(True, wavelengths=[2.], emission_max_escape_optical_depth=5.)
    m2._read_monochromatic(v)
    assert m2._monochromatic_emission_max_escape_optical_depth == value
