
from copy import deepcopy

import numpy as np
from numpy.testing import assert_equal
import pytest
from astropy import units as u

from .. import Model, AnalyticalYSOModel, ModelOutput
from ...grid import AMRGrid
from ...util.functions import virtual_file
from .test_helpers import random_id, get_test_dust, assert_identical_results


def test_io_minimal(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = Model()
    m1.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m1.set_n_initial_iterations(0)
    m1.set_n_photons(imaging=10)
    m1.write(filename)

    m2 = Model.read(filename)


def test_io_minimal_mono(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = Model()
    m1.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m1.set_n_initial_iterations(0)
    m1.set_monochromatic(True, wavelengths=[1, 4.5, 7.7])
    m1.set_n_photons(imaging_sources=10, imaging_dust=10)
    m1.write(filename)

    m2 = Model.read(filename)

    assert m1._monochromatic == m2._monochromatic
    assert_equal(m1._frequencies, m2._frequencies)


@pytest.mark.parametrize(('value'), [True, False])
def test_io_monochromatic(value):
    m1 = Model()
    if value:
        m1.set_monochromatic(value, wavelengths=[1, 4.5, 7.7])
    else:
        m1.set_monochromatic(value)
    v = virtual_file()
    m1._write_monochromatic(v)
    m2 = Model()
    m2._read_monochromatic(v)
    assert m1._monochromatic == m2._monochromatic
    if value:
        assert_equal(m1._frequencies, m2._frequencies)


def test_io_monochromatic_full(tmpdir):
    """
    Regression for issues described in #77 - namely that when reading and
    writing out a model with monochromatic settings, images/SEDs would not work
    properly (raising an exception).
    """

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    i = m.add_peeled_images(image=True)
    i.set_viewing_angles([45.], [45.])
    i.set_image_size(256, 256)
    i.set_image_limits(-1., 1., -1., 1.)

    m.set_monochromatic(True, wavelengths=[0.1, 1., 10.])

    m.set_n_photons(initial=1000, imaging_sources=1000, imaging_dust=1000)

    m.write(filename1)

    m2 = Model.read(filename1)
    m2.write(filename2)


def test_io_monochromatic_full_2(tmpdir):
    """
    Regression for issues described in #77 - namely that when reading and
    writing out a model with monochromatic settings, images/SEDs would not work
    properly (raising an exception).

    Changing the order to setting before the monochromatic mode is set changes
    how things are handled internally, hence the similar test.
    """

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    m.set_monochromatic(True, wavelengths=[0.1, 1., 10.])

    i = m.add_peeled_images(image=True)
    i.set_viewing_angles([45.], [45.])
    i.set_image_size(256, 256)
    i.set_image_limits(-1., 1., -1., 1.)

    m.set_n_photons(initial=1000, imaging_sources=1000, imaging_dust=1000)

    m.write(filename1)

    m2 = Model.read(filename1)
    m2.write(filename2)


def test_io_voronoi(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    x = [1., 3., 2., 5., 2.]
    y = [3., 2., 8., 3., 6.]
    z = [4., 3., 1., 1., 2.]

    m1 = Model()
    m1.set_voronoi_grid(x, y, z)
    m1.set_n_initial_iterations(0)
    m1.set_n_photons(imaging=10)
    m1.write(filename)

    m2 = Model.read(filename)


@pytest.mark.requires_hyperion_binaries
def test_read_then_monochromatic(tmpdir):

    # If we run a model then read in the output and make a new model with
    # monochromatic output, Hyperion used to raise an error that the number
    # of non-monochromatic had already been set by there was no way to
    # fix that with the public API.

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    s = m.add_point_source()
    s.position = (0.4, 0., 0.)
    s.luminosity = 1000
    s.temperature = 6000.

    m.set_n_photons(initial=10, imaging=0)

    m.write(filename1 + '.rtin')
    m.run(filename1 + '.rtout')

    m2 = Model.read(filename1 + '.rtout')

    m.set_monochromatic(True, wavelengths=[1])
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging_sources=10, imaging_dust=10)

    m.write(filename2 + '.rtin', overwrite=True)


ALL_GRID_TYPES = ['car', 'cyl', 'sph', 'amr', 'oct', 'vor']


def set_test_grid(model, grid_type):
    if grid_type == 'car':
        model.set_cartesian_grid([-1., 0., 1.], [-2., 2.], [-3., -1., 3.])
    elif grid_type == 'cyl':
        model.set_cylindrical_polar_grid([0., 0.5, 1.], [-1., 1.], [0., np.pi, 2. * np.pi])
    elif grid_type == 'sph':
        model.set_spherical_polar_grid([0., 0.5, 1.], [0., np.pi / 2., np.pi], [0., 2. * np.pi])
    elif grid_type == 'amr':
        amr = AMRGrid()
        level = amr.add_level()
        grid = level.add_grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 4, 4, 4
        model.set_amr_grid(amr)
    elif grid_type == 'oct':
        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        model.set_octree_grid(0., 0., 0., 10., 10., 10., np.array(refined).astype(bool))
    elif grid_type == 'vor':
        model.set_voronoi_grid([1., 3., 2., 5., 2.],
                               [3., 2., 8., 3., 6.],
                               [4., 3., 1., 1., 2.])


def random_quantity(grid, random_state):
    # Return an array (or AMR grid view) of random values with the shape of
    # the grid
    if isinstance(grid, AMRGrid):
        amr = deepcopy(grid)
        for level in amr.levels:
            for g in level.grids:
                g.quantities['tmp'] = random_state.uniform(0.5, 1., (g.nz, g.ny, g.nx))
        return amr['tmp']
    else:
        return random_state.uniform(0.5, 1., grid.shape)


def get_full_model(grid_type, specific_energy=True, map_source=True, output_config=True):

    random_state = np.random.RandomState(12345)

    m = Model()

    set_test_grid(m, grid_type)

    # Two dust populations, with density and specific energy

    for idust in range(2):
        m.add_density_grid(random_quantity(m.grid, random_state), get_test_dust(),
                           specific_energy=random_quantity(m.grid, random_state) if specific_energy else None)

    m.set_minimum_specific_energy([0.1, 0.2])

    # Sources of all types

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 5000.
    s.position = (0.1, 0.2, 0.3)

    s = m.add_point_source_collection()
    s.luminosity = np.array([1., 2.])
    s.position = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
    s.spectrum = ([1.e12, 1.e13, 1.e14], [1., 2., 3.])

    s = m.add_spherical_source()
    s.luminosity = 2.
    s.radius = 0.1
    s.temperature = 6000.
    s.limb = True
    spot = s.add_spot()
    spot.luminosity = 0.1
    spot.temperature = 7000.
    spot.longitude = 20.
    spot.latitude = 30.
    spot.radius = 10.

    s = m.add_external_spherical_source()
    s.luminosity = 3.
    s.radius = 100.
    s.temperature = 1000.

    s = m.add_external_box_source()
    s.luminosity = 4.
    s.bounds = [[-100., 100.], [-200., 200.], [-300., 300.]]
    s.temperature = 2000.

    if map_source:
        s = m.add_map_source()
        s.luminosity = 5.
        s.map = random_quantity(m.grid, random_state)
        s.temperature = 3000.

    s = m.add_plane_parallel_source()
    s.luminosity = 6.
    s.radius = 0.5
    s.position = (0.1, 0.2, 0.3)
    s.direction = (30., 40.)
    s.temperature = 4000.

    # Peeled and binned images

    i = m.add_peeled_images(sed=True, image=True)
    i.set_viewing_angles([10., 20.], [30., 40.])
    i.set_image_size(4, 6)
    i.set_image_limits(-1., 1., -2., 2.)
    i.set_aperture_radii(3, 0.1, 1.)
    i.set_wavelength_range(5, 0.1, 100.)
    i.set_track_origin('detailed')
    i.set_uncertainties(True)
    i.set_stokes(True)
    i.set_peeloff_origin((0.1, 0.2, 0.3))
    i.set_output_bytes(4)

    i = m.add_peeled_images(sed=False, image=True)
    i.set_inside_observer((0.1, 0.2, 0.3))
    i.set_viewing_angles([10.], [20.])
    i.set_image_size(4, 4)
    i.set_image_limits(10., -10., -10., 10.)
    i.set_wavelength_range(3, 1., 10.)
    i.set_ignore_optical_depth(True)
    i.set_depth(0.1, 0.5)

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([45.], [45.])
    i.set_track_origin('scatterings', n_scat=3)
    f = i.add_filter()
    f.name = 'F1'
    f.spectral_coord = [1, 1.1, 1.2, 1.3] * u.micron
    f.transmission = [0., 100., 50, 0.] * u.percent
    f.detector_type = 'photons'
    f.alpha = 0.
    f.central_spectral_coord = 1.15 * u.micron

    i = m.add_binned_images(sed=True, image=True)
    i.set_viewing_bins(3, 4)
    i.set_image_size(5, 5)
    i.set_image_limits(-1., 1., -1., 1.)
    i.set_wavelength_range(4, 0.1, 10.)

    # Non-default run configuration

    m.set_seed(-12345)
    m.set_propagation_check_frequency(0.1)
    m.set_n_initial_iterations(3)
    m.set_raytracing(True)
    m.set_n_photons(initial=100, imaging=200, raytracing_sources=10,
                    raytracing_dust=20, stats=50)
    m.set_max_interactions(1000, warn=False)
    m.set_max_reabsorptions(2000, warn=False)
    m.set_max_peeloff_optical_depth(10.)
    # PDA and MRW are not implemented for Voronoi grids, and writing such a
    # model raises an error
    if grid_type != 'vor':
        m.set_pda(True)
        m.set_mrw(True, gamma=2., inter_max=500, warn=False)
    m.set_convergence(True, percentile=95., absolute=2., relative=1.1)
    m.set_kill_on_absorb(True)
    m.set_forced_first_interaction(False)
    m.set_compute_escape_column_density(True, rays_per_cell=5)
    m.set_output_bytes(4)
    m.set_sample_sources_evenly(True)
    m.set_enforce_energy_range(False)
    m.set_copy_input(False)
    if specific_energy:
        m.set_specific_energy_type('additional')

    # Non-default output configuration

    if output_config:
        m.conf.output.output_density = 'all'
        m.conf.output.output_density_diff = 'last'
        m.conf.output.output_specific_energy = 'all'
        m.conf.output.output_n_photons = 'last'

    return m


def assert_attributes_equal(obj1, obj2, exclude=()):
    # Check that all the attributes of obj1 are present and equal in obj2.
    # The object that was read in can have extra attributes that are not
    # relevant (e.g. aperture settings for image-only configurations).
    keys1 = set(vars(obj1)) - set(exclude) - {'_attributes', '_frozen', '_final'}
    assert keys1 <= set(vars(obj2))
    for key in sorted(keys1):
        value1, value2 = getattr(obj1, key), getattr(obj2, key)
        if isinstance(value1, (list, tuple)) and len(value1) > 0 and hasattr(value1[0], '__dict__'):
            assert len(value1) == len(value2), key
            for item1, item2 in zip(value1, value2):
                assert_attributes_equal(item1, item2)
        elif hasattr(value1, '__dict__') and not isinstance(value1, np.ndarray):
            assert_attributes_equal(value1, value2)
        else:
            assert_equal(value1, value2, err_msg=key)


def get_quantity_arrays(grid, quantity):
    # Return a list of the arrays for each dust population, with the arrays
    # for all AMR grids concatenated
    if isinstance(grid, AMRGrid):
        arrays = [g.quantities[quantity] for level in grid.levels for g in level.grids]
        return [np.hstack([np.ravel(a[idust]) for a in arrays]) for idust in range(len(arrays[0]))]
    else:
        return [np.asarray(array) for array in grid.quantities[quantity]]


GRID_ATTRIBUTES = {'car': ('x_wall', 'y_wall', 'z_wall'),
                   'cyl': ('w_wall', 'z_wall', 'p_wall'),
                   'sph': ('r_wall', 't_wall', 'p_wall'),
                   'oct': ('x', 'y', 'z', 'dx', 'dy', 'dz', 'refined'),
                   'vor': ('_x', '_y', '_z', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')}

AMR_GRID_ATTRIBUTES = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'nx', 'ny', 'nz')


@pytest.mark.parametrize('grid_type', ALL_GRID_TYPES)
def test_io_full_roundtrip(tmpdir, grid_type):

    # Write out a model with all features set, read it back in, and check
    # that the model is the same, and that writing it out again gives the
    # same file. The features with known round-tripping bugs are tested
    # separately below.

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m1 = get_full_model(grid_type,
                        specific_energy=grid_type != 'amr',
                        map_source=grid_type != 'amr',
                        output_config=False)
    m1.write(filename1)

    m2 = Model.read(filename1)

    # Geometry
    assert type(m2.grid) is type(m1.grid)
    if grid_type == 'amr':
        assert len(m2.grid.levels) == len(m1.grid.levels)
        for level1, level2 in zip(m1.grid.levels, m2.grid.levels):
            assert len(level2.grids) == len(level1.grids)
            for grid1, grid2 in zip(level1.grids, level2.grids):
                for attr in AMR_GRID_ATTRIBUTES:
                    assert getattr(grid2, attr) == getattr(grid1, attr)
    else:
        assert m2.grid.shape == m1.grid.shape
        for attr in GRID_ATTRIBUTES[grid_type]:
            assert_equal(getattr(m2.grid, attr), getattr(m1.grid, attr), err_msg=attr)

    # Quantities
    for quantity in ('density', 'specific_energy'):
        assert (quantity in m2.grid) is (quantity in m1.grid)
        if quantity in m1.grid:
            arrays1 = get_quantity_arrays(m1.grid, quantity)
            arrays2 = get_quantity_arrays(m2.grid, quantity)
            assert len(arrays2) == len(arrays1) == 2
            for array1, array2 in zip(arrays1, arrays2):
                assert_equal(array2, array1)
    assert_equal(m2._minimum_specific_energy, m1._minimum_specific_energy)

    # Dust
    assert len(m2.dust) == len(m1.dust)

    # Sources
    assert [type(s) for s in m2.sources] == [type(s) for s in m1.sources]
    for s1, s2 in zip(m1.sources, m2.sources):
        assert_attributes_equal(s1, s2)

    # Run and output configuration
    assert_attributes_equal(m1.conf, m2.conf)

    # Image configuration
    assert_attributes_equal(m1.binned_output, m2.binned_output)
    assert len(m2.peeled_output) == len(m1.peeled_output)
    for i1, i2 in zip(m1.peeled_output, m2.peeled_output):
        # n_view is tested separately below. Filters are normalized and
        # converted to frequency when written, so we only check the scalar
        # attributes here.
        assert_attributes_equal(i1, i2, exclude=('n_view', '_filters'))
        assert len(i2._filters) == len(i1._filters)
        for f1, f2 in zip(i1._filters, i2._filters):
            assert f2.name == f1.name
            assert f2.alpha == f1.alpha
            assert f2.detector_type == f1.detector_type

    # Writing out the model again should give identical HDF5 content
    m2.write(filename2)
    assert_identical_results(filename1, filename2)


@pytest.mark.xfail(strict=True, reason="Model.use_output_config discards the result of OutputConf.read")
def test_io_output_config(tmpdir):

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m1 = get_full_model('car')
    m1.write(filename1)

    m2 = Model.read(filename1)
    assert_attributes_equal(m1.conf.output, m2.conf.output)

    m2.write(filename2)
    assert_identical_results(filename1, filename2)


@pytest.mark.xfail(strict=True, reason="PeeledImageConf._read_viewing_angles does not set n_view")
def test_io_peeled_n_view(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = get_full_model('car')
    m1.write(filename)

    m2 = Model.read(filename)
    assert [i.n_view for i in m2.peeled_output] == [i.n_view for i in m1.peeled_output]


@pytest.mark.xfail(strict=True, reason="Model.use_quantities only looks for specific_energy at the top level of the Quantities group, so misses it for AMR grids")
def test_io_amr_specific_energy(tmpdir):

    filename = tmpdir.join(random_id()).strpath

    m1 = get_full_model('amr', map_source=False)
    m1.write(filename)

    m2 = Model.read(filename)
    assert 'specific_energy' in m2.grid
    for array1, array2 in zip(get_quantity_arrays(m1.grid, 'specific_energy'),
                              get_quantity_arrays(m2.grid, 'specific_energy')):
        assert_equal(array2, array1)


@pytest.mark.xfail(strict=True, raises=KeyError, reason="MapSource.read expects a single dataset, but AMR maps are written as a group of levels")
def test_io_amr_map_source(tmpdir):

    filename1 = tmpdir.join(random_id()).strpath
    filename2 = tmpdir.join(random_id()).strpath

    m1 = get_full_model('amr', specific_energy=False, output_config=False)
    m1.write(filename1)

    m2 = Model.read(filename1)
    m2.write(filename2)
    assert_identical_results(filename1, filename2)


def test_analytical_yso_model_read():
    with pytest.raises(Exception):
        AnalyticalYSOModel.read('model.rtin')


@pytest.mark.requires_hyperion_binaries
def test_io_read_only_initial(tmpdir):

    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    m = Model()
    m.set_cartesian_grid([-1., 0., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1., 2.]]]), get_test_dust())
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 5000.
    m.set_n_photons(initial=1000, imaging=0)
    m.set_n_initial_iterations(2)
    m.write(input_file)
    m.run(output_file)

    # By default, only the input quantities are read
    m_initial = Model.read(output_file)
    assert 'specific_energy' not in m_initial.grid
    assert_equal(m_initial.grid['density'][0].array, m.grid['density'][0].array)

    # Otherwise, the specific energy from the last iteration should be read
    m_final = Model.read(output_file, only_initial=False)
    expected = ModelOutput(output_file).get_quantities(iteration=-1)
    assert_equal(m_final.grid['specific_energy'][0].array,
                 expected['specific_energy'][0].array)
    assert np.all(m_final.grid['specific_energy'][0].array > 0.)
