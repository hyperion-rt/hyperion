# The purpose of these tests is to check that the results produced by
# Hyperion agree at the bit-level. The data for these tests will need to be
# changed if algorithms are updated, and bugs are fixed, but otherwise there
# is no reason we should expect a change in results from one commit to the
# next. Since these files take place, we should minimize the number of tests
# to run.

import os
import shutil
import itertools

from astropy.tests.helper import pytest
import numpy as np

from .test_helpers import random_id, assert_identical_results
from .. import Model, AnalyticalYSOModel
from ...util.constants import pc, lsun, c, au, msun, pi, sigma, rsun
from ...grid import CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid, AMRGrid, OctreeGrid
from ...dust import IsotropicDust, SphericalDust

GRID_TYPES = ['car', 'cyl', 'sph', 'amr', 'oct']

DATA = os.path.join(os.path.dirname(__file__), 'data')

bit_level = pytest.mark.skipif(str(not pytest.config.getoption('enable_bit_level_tests')))


@pytest.fixture(scope="module")
def generate(request):
    generate_reference = request.config.getvalue("generate_reference")
    if generate_reference is None:
        return False
    else:
        return generate_reference


def setup_all_grid_types(self, u, d):
    '''
    All grids are guaranteed to cover the volume from -u to u in x, y, z
    '''

    np.random.seed(141412)

    self.grid = {}

    # Cartesian
    x = np.linspace(-u, u, 8)
    y = np.linspace(-u, u, 6)
    z = np.linspace(-u, u, 4)
    self.grid['car'] = CartesianGrid(x, y, z)

    # Cylindrical polar
    w = np.linspace(0., 2. * u, 8)
    z = np.linspace(-u, u, 4)
    p = np.linspace(0., 2. * np.pi, 6)
    self.grid['cyl'] = CylindricalPolarGrid(w, z, p)

    # Spherical polar
    r = np.linspace(0., 3. * u, 6)
    t = np.linspace(0., np.pi, 8)
    p = np.linspace(0., 2. * np.pi, 4)
    self.grid['sph'] = SphericalPolarGrid(r, t, p)

    # AMR

    self.grid['amr'] = AMRGrid()

    level1 = self.grid['amr'].add_level()
    grid = level1.add_grid()
    grid.xmin, grid.xmax = -u, u
    grid.ymin, grid.ymax = -u, u
    grid.zmin, grid.zmax = -u, u
    grid.nx, grid.ny, grid.nz = 8, 6, 4
    grid.quantities['density'] = np.random.random((4, 6, 8)) * d
    grid.quantities['density_2'] = np.random.random((4, 6, 8)) * d
    grid.quantities['density_3'] = np.random.random((4, 6, 8)) * d

    level2 = self.grid['amr'].add_level()
    grid2 = level2.add_grid()
    grid2.xmin, grid2.xmax = -u, 0.
    grid2.ymin, grid2.ymax = -u, 0.
    grid2.zmin, grid2.zmax = -u, 0.
    grid2.nx, grid2.ny, grid2.nz = 4, 6, 20
    grid2.quantities['density'] = np.random.random((20, 6, 4)) * d
    grid2.quantities['density_2'] = np.random.random((20, 6, 4)) * d
    grid2.quantities['density_3'] = np.random.random((20, 6, 4)) * d

    # Octree
    refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    self.grid['oct'] = OctreeGrid(0., 0., 0., u, u, u, np.array(refined).astype(bool))

    # Set up initial densities
    self.density = {}
    self.density['car'] = np.random.random(self.grid['car'].shape) * d
    self.density['cyl'] = np.random.random(self.grid['cyl'].shape) * d
    self.density['sph'] = np.random.random(self.grid['sph'].shape) * d
    self.density['amr'] = self.grid['amr']['density']
    self.density['oct'] = np.random.random(len(refined)) * d

    # Second set of densities
    self.density_2 = {}
    self.density_2['car'] = np.random.random(self.grid['car'].shape) * d
    self.density_2['cyl'] = np.random.random(self.grid['cyl'].shape) * d
    self.density_2['sph'] = np.random.random(self.grid['sph'].shape) * d
    self.density_2['amr'] = self.grid['amr']['density_2']
    self.density_2['oct'] = np.random.random(len(refined)) * d

    # Third set of densities
    self.density_3 = {}
    self.density_3['car'] = np.random.random(self.grid['car'].shape) * d
    self.density_3['cyl'] = np.random.random(self.grid['cyl'].shape) * d
    self.density_3['sph'] = np.random.random(self.grid['sph'].shape) * d
    self.density_3['amr'] = self.grid['amr']['density_3']
    self.density_3['oct'] = np.random.random(len(refined)) * d


def function_name():
    import sys
    import inspect
    caller = sys._getframe(1)
    args, _, _, values = inspect.getargvalues(caller)
    name = [caller.f_code.co_name]
    for arg in args:
        if arg not in ['self', 'generate', 'tmpdir']:
            name += ["{0}={1}".format(arg, values[arg])]
    name = '.'.join(name)
    return name


class TestBasic(object):

    def setup_class(self):
        setup_all_grid_types(self, pc, 1.e-20)
        self.dust_file = os.path.join(DATA, 'kmh_lite.hdf5')

    @bit_level
    @pytest.mark.parametrize(('grid_type', 'sample_sources_evenly', 'multiple_densities'), list(itertools.product(GRID_TYPES, [False, True], [False, True])))
    def test_specific_energy(self, tmpdir, grid_type, sample_sources_evenly, multiple_densities, generate):

        np.random.seed(12345)

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust_file)
        if multiple_densities:
            m.add_density_grid(self.density_2[grid_type], self.dust_file)
            m.add_density_grid(self.density_3[grid_type], self.dust_file)

        for i in range(5):
            s = m.add_point_source()
            s.luminosity = np.random.random() * lsun
            s.temperature = np.random.uniform(2000., 10000.)
            s.position = np.random.uniform(-pc, pc, 3)

        m.set_n_photons(initial=10000, imaging=0)

        m.set_sample_sources_evenly(sample_sources_evenly)

        m.conf.output.output_specific_energy = 'all'

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)

    @bit_level
    @pytest.mark.parametrize(('grid_type', 'raytracing', 'sample_sources_evenly'), list(itertools.product(GRID_TYPES, [False, True], [False, True])))
    def test_peeloff(self, tmpdir, grid_type, raytracing, sample_sources_evenly, generate):

        np.random.seed(12345)

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust_file)

        for i in range(5):
            s = m.add_point_source()
            s.luminosity = np.random.random() * lsun
            s.temperature = np.random.uniform(2000., 10000.)
            s.position = np.random.uniform(-pc, pc, 3)

        m.set_raytracing(raytracing)
        if raytracing:
            m.set_n_photons(initial=1000, imaging=5000, raytracing_sources=2000, raytracing_dust=3000)
        else:
            m.set_n_photons(initial=1000, imaging=5000)

        m.set_sample_sources_evenly(sample_sources_evenly)

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(5, 0.05, 200.)
        i_p.set_viewing_angles([33.4, 110.], [65.4, 103.2])
        i_p.set_image_size(4, 5)
        i_p.set_image_limits(-0.8 * pc, 0.8 * pc, -pc, pc)
        i_p.set_aperture_range(5, 0.1 * pc, pc)
        i_p.set_stokes(True)

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(4, 0.05, 200.)
        i_p.set_viewing_angles([22.1], [203.2])
        i_p.set_image_size(6, 6)
        i_p.set_image_limits(-pc, pc, -pc, pc)
        i_p.set_aperture_range(2, 0.5 * pc, pc)
        i_p.set_track_origin('basic')
        i_p.set_stokes(True)

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(4, 0.05, 200.)
        i_p.set_viewing_angles([22.1], [203.2])
        i_p.set_image_size(6, 6)
        i_p.set_image_limits(-pc, pc, -pc, pc)
        i_p.set_aperture_range(2, 0.5 * pc, pc)
        i_p.set_track_origin('detailed')
        i_p.set_stokes(True)

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)


class TestPascucciBenchmark(object):
    '''
    Run a very low signal-to-noise version of the Pascucci benchmark models.
    The point is not to have high signal-to-noise but to look for bit-level
    changes.
    '''

    def setup_class(self):

        optSi = '''
           0.1200000       5.8811883E-14   1.1439794E-13
           0.1400000       5.8397304E-14   1.1660481E-13
           0.1600000       6.2787212E-14   1.2265337E-13
           0.1800000       5.3791878E-14   1.1174947E-13
           0.2000000       6.5517043E-14   1.0667109E-13
           0.2150000       1.0607825E-13   1.4219348E-13
           0.2200000       1.0908588E-13   1.3783945E-13
           0.2300000       1.3212733E-13   1.5250003E-13
           0.2500000       1.8018174E-13   2.0492832E-13
           0.2740000       1.6293549E-13   1.8044428E-13
           0.3000000       1.8149981E-13   1.9903920E-13
           0.3440000       1.4926875E-13   1.6465974E-13
           0.4000000       1.2503861E-13   1.3638487E-13
           0.4400000       1.1388763E-13   1.2583496E-13
           0.5500000       5.3850835E-14   6.1417044E-14
           0.7000000       2.4657287E-14   2.9103080E-14
           0.9000000       9.7663111E-15   1.2698527E-14
            1.100000       4.3912416E-15   6.5586447E-15
            1.400000       1.6462753E-15   3.2172486E-15
            1.650000       8.4103368E-16   2.1375954E-15
            2.000000       3.8364113E-16   1.4456605E-15
            2.200000       2.6030839E-16   1.2321842E-15
            2.600000       1.3101047E-16   9.7629208E-16
            3.000000       7.2728141E-17   8.3056811E-16
            3.200000       5.5601797E-17   7.8231688E-16
            3.600000       3.3966690E-17   7.0568777E-16
            4.000000       2.1781298E-17   6.4674454E-16
            5.000000       8.0873324E-18   5.5537600E-16
            6.000000       3.2988649E-18   5.4991692E-16
            6.280000       2.5180852E-18   5.5436766E-16
            6.300000       2.4702033E-18   5.5468302E-16
            6.320000       2.4232568E-18   5.5500172E-16
            6.500000       2.0396911E-18   5.5833034E-16
            8.000000       2.8743327E-19   1.7708248E-15
            9.500000       1.5971045E-18   7.2037611E-15
            10.00000       1.6023445E-18   6.5482203E-15
            11.50000       9.9434658E-19   3.7331190E-15
            11.51500       9.8820641E-19   3.7104210E-15
            11.52500       9.8415064E-19   3.6953642E-15
            11.54000       9.7812174E-19   3.6728894E-15
            12.00000       8.2106490E-19   3.0439779E-15
            14.00000       3.1000245E-19   1.5214940E-15
            16.00000       1.5696574E-19   2.0153394E-15
            18.00000       1.4169530E-19   2.5528336E-15
            20.00000       1.2279225E-19   2.2793345E-15
            24.00000       6.5145023E-20   1.5326092E-15
            27.50000       3.9240954E-20   1.1453648E-15
            32.50000       2.0514811E-20   8.2030043E-16
            37.50000       1.1663588E-20   6.2319899E-16
            45.00000       5.7365085E-21   4.2517091E-16
            55.00000       2.6121877E-21   2.7811657E-16
            70.00000       1.0024304E-21   1.6293315E-16
            90.00000       3.6573556E-22   9.8668075E-17
            110.0000       1.6381540E-22   6.4738834E-17
            135.0000       7.2083289E-23   4.2774677E-17
            175.0000       2.5504432E-23   2.4758561E-17
            250.0000       6.1052632E-24   9.8178232E-18
            400.0000       9.3257302E-25   4.8822549E-18
            700.0000       9.9398631E-26   1.4894551E-18
            1200.000       1.1479267E-26   5.0961962E-19
            2000.000       1.4912529E-27   1.8238127E-19
            '''

        # Read in dust
        from StringIO import StringIO
        data = np.loadtxt(StringIO(optSi), comments=';',
                          dtype=[('wav', float), ('csca', float),
                                 ('cext', float)])

        # Convert SI to CGS
        data['cext'] *= 1.e4
        data['csca'] *= 1.e4
        grain_size = 0.12 * 1.e-4
        grain_density = 3.6
        wav = data['wav']
        nu = c / (wav * 1.e-4)
        chi = data['cext'] / (4. * pi / 3. * grain_size ** 3. * grain_density)
        albedo = data['csca'] / data['cext']

        # Set up dust object
        import tempfile
        self.tmpdir = tempfile.mkdtemp()
        self.dust_file = os.path.join(self.tmpdir, random_id())
        dust = IsotropicDust(nu[::-1], albedo[::-1], chi[::-1])
        dust.optical_properties.extrapolate_wav(1.e-3, 1.e5)
        dust.set_lte_emissivities(n_temp=100, temp_min=0.1, temp_max=1600.)
        dust.write(self.dust_file)

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    @bit_level
    @pytest.mark.parametrize(('tau'), [0.1, 1, 10, 100])
    def test_pascucci(self, tmpdir, tau, generate):

        print(generate)

        m = AnalyticalYSOModel()

        # Star:
        #
        # T_star = 5,800K
        # R_star = R_sun
        # M_star = M_sun

        m.star.radius = 1.  # = 1 cm (i.e. point source)
        m.star.temperature = 5800.
        m.star.luminosity = 4. * pi * rsun ** 2 * sigma * 5800. ** 4

        # Disk:
        #
        # rho(r,z) = rho_0 * (r/rd)^-1.0 * exp(-pi * (z/(2*h(r)))^2)
        # h(r) = zd * (r/rd)^1.125
        # rd = rmax / 2.
        # zd = rmax / 8.
        # r_in = 1AU, r_out = 1000AU\
        # Sharp vertical edges in cylindrical polar coordinates
        # Disk mass (dust!) = 1.1e-7, 1.1e-6, 1.1e-5, 1.1e-4

        disk = m.add_flared_disk()
        disk.p = 0.125
        disk.beta = 1.125
        disk.mass = 1.113838e-6 * msun * tau
        disk.rmin = 1. * au
        disk.rmax = 1000. * au
        disk.h_0 = 125 * au * np.sqrt(2. / pi)
        disk.r_0 = 500 * au

        # Dust:
        #
        # Spherical grains, 1 micron, 3.5g/cm^3

        disk.dust = self.dust_file

        # SEDs/Images:
        #
        # cos(i) = 0.05 to 0.95 in steps of 0.1
        # Images computed at 1 micron
        # 251 pixels/900 AU across, at 140 pc

        theta = np.array([12.5, 42.5, 77.5])
        phi = np.array([30.0, 30.0, 30.0])

        image = m.add_peeled_images()
        image.set_viewing_angles(theta, phi)
        image.set_image_size(1, 1)
        image.set_image_limits(-1500. * au, 1500. * au, -1500. * au, 1500. * au)
        image.set_aperture_range(1, 1500. * au, 1500. * au)
        image.set_wavelength_range(61, 1, 61)
        image.set_stokes(True)

        m.set_raytracing(True)
        m.set_n_initial_iterations(5)
        # Don't test for convergence, since this is a low signal-to-noise models

        # Use a lower-resolution grid
        m.set_spherical_polar_grid_auto(100, 30, 1, rmax=1300. * au)

        wavelengths = [0.12, 0.14, 0.16, 0.18, 0.2, 0.215, 0.22, 0.23, 0.25,
                       0.274, 0.3, 0.344, 0.4, 0.44, 0.55, 0.7, 0.9, 1.1,
                       1.4, 1.65, 2, 2.2, 2.6, 3, 3.2, 3.6, 4, 5, 6, 6.28,
                       6.3, 6.32, 6.5, 8, 9.5, 10, 11.5, 11.515016,
                       11.524977, 11.540016, 12, 14, 16, 18, 20, 24, 27.5,
                       32.5, 37.5, 45, 55, 70, 90, 110, 135, 175, 250, 400,
                       700, 1200, 2000]

        m.set_monochromatic(True, wavelengths=wavelengths)

        m.set_n_photons(initial=1000, imaging_sources=1000, imaging_dust=1000,
                        raytracing_sources=1000, raytracing_dust=1000)

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)


class TestPinteBenchmark(object):
    '''
    Run a very low signal-to-noise version of the Pinte benchmark models.
    The point is not to have high signal-to-noise but to look for bit-level
    changes.

    Notes
    -----

    The current tests do not test the imaging part of the Pinte benchmark.
    '''

    @bit_level
    @pytest.mark.parametrize(('tau'), [1000, 10000, 100000, 1000000])
    def test_pinte_seds(self, tmpdir, tau, generate):

        m = AnalyticalYSOModel()

        # Star:
        #
        # T_star = 4,000K
        # R_star = 2 * R_sun

        m.star.radius = 2. * rsun
        m.star.temperature = 4000.
        m.star.luminosity = 4. * pi * (2. * rsun) ** 2. * sigma * 4000. ** 4.

        # Disk:
        #
        # Sigma(r) = Sigma_0 * (r/r_0)^-1.5
        # h(r) = h_0 * (r/r_0)^1.125
        # h_0 = 10 AU at 100AU
        # r_in = 0.1AU, r_out = 400AU\
        # Sharp vertical edges in cylindrical polar coordinates
        # Disk mass (dust!) = 3e-8, 3e-7, 3e-6, 3e-5

        disk = m.add_flared_disk()
        disk.p = -1.5
        disk.beta = 1.125
        disk.mass = 3.e-8 * msun * tau / 1.e3
        disk.rmin = 0.1 * au
        disk.rmax = 400 * au
        disk.h_0 = 10 * au
        disk.r_0 = 100. * au

        disk.cylindrical_inner_rim = True
        disk.cylindrical_outer_rim = True

        # Dust:
        #
        # Spherical grains, 1 micron, 3.5g/cm^3

        disk.dust = SphericalDust(os.path.join(DATA, 'pinte_dust_lite.hdf5'))

        # SEDs/Images:
        #
        # cos(i) = 0.05 to 0.95 in steps of 0.1
        # Images computed at 1 micron
        # 251 pixels/900 AU across, at 140 pc

        # theta = np.degrees(np.arccos(np.linspace(0.05,0.95,10)))
        # phi = np.ones(10) * 30

        theta = np.degrees(np.arccos(np.array([0.95, 0.25, 0.15, 0.05])))
        phi = np.array([45., 45., 45., 45.])

        image = m.add_peeled_images()
        image.set_viewing_angles(theta, phi)
        image.set_image_size(1, 1)
        image.set_image_limits(-450. * au, 450. * au, -450. * au, 450. * au)
        image.set_aperture_range(1, 450. * au, 450. * au)
        image.set_wavelength_range(2000, 0.01, 5000.)
        image.set_stokes(True)

        m.set_raytracing(True)

        m.set_n_initial_iterations(10)
        m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

        m.set_cylindrical_polar_grid_auto(100, 30, 1)

        wavelengths = [0.110635, 0.135419, 0.165755, 0.202887, 0.248336,
                       0.303967, 0.372060, 0.455408, 0.557426, 0.682297,
                       0.835142, 1.02223, 1.25122, 1.53151, 1.87459,
                       2.29453, 2.80854, 3.43769, 4.20779, 5.15039, 6.30416,
                       7.71638, 9.44497, 11.5608, 14.1506, 17.3205, 21.2006,
                       25.9498, 31.7629, 38.8783, 47.5876, 58.2480, 71.2964,
                       87.2678, 106.817, 130.746, 160.035, 195.885, 239.766,
                       293.477, 359.220, 439.691, 538.188, 658.751, 806.321,
                       986.948, 1208.04, 1478.66, 1809.90, 2215.34, 2711.61]

        m.set_monochromatic(True, wavelengths=wavelengths)

        m.set_mrw(True, gamma=2.)
        # Don't use the PDA here because it's too slow when there are too few photons

        m.set_n_photons(initial=5000, imaging_sources=100, imaging_dust=2000,
                        raytracing_sources=1000, raytracing_dust=1000)

        m.set_max_interactions(10000000)

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)

    @bit_level
    @pytest.mark.parametrize(('tau'), [1000, 10000, 100000, 1000000])
    def test_pinte_images(self, tmpdir, tau, generate):

        m = AnalyticalYSOModel()

        # Star:
        #
        # T_star = 4,000K
        # R_star = 2 * R_sun

        m.star.radius = 2. * rsun
        m.star.temperature = 4000.
        m.star.luminosity = 4. * pi * (2. * rsun) ** 2. * sigma * 4000. ** 4.

        # Disk:
        #
        # Sigma(r) = Sigma_0 * (r/r_0)^-1.5
        # h(r) = h_0 * (r/r_0)^1.125
        # h_0 = 10 AU at 100AU
        # r_in = 0.1AU, r_out = 400AU\
        # Sharp vertical edges in cylindrical polar coordinates
        # Disk mass (dust!) = 3e-8, 3e-7, 3e-6, 3e-5

        disk = m.add_flared_disk()
        disk.p = -1.5
        disk.beta = 1.125
        disk.mass = 3.e-8 * msun * tau / 1.e3
        disk.rmin = 0.1 * au
        disk.rmax = 400 * au
        disk.h_0 = 10 * au
        disk.r_0 = 100. * au

        disk.cylindrical_inner_rim = True
        disk.cylindrical_outer_rim = True

        # Dust:
        #
        # Spherical grains, 1 micron, 3.5g/cm^3

        disk.dust = SphericalDust(os.path.join(DATA, 'pinte_dust_lite.hdf5'))

        # SEDs/Images:
        #
        # cos(i) = 0.05 to 0.95 in steps of 0.1
        # Images computed at 1 micron
        # 251 pixels/900 AU across, at 140 pc

        # theta = np.degrees(np.arccos(np.linspace(0.05,0.95,10)))
        # phi = np.ones(10) * 30

        theta = np.array([69.5, 87.1])
        phi = np.array([45., 45.])

        image = m.add_peeled_images()
        image.set_viewing_angles(theta, phi)
        image.set_image_size(51, 51)
        image.set_image_limits(-450. * au, 450. * au, -450. * au, 450. * au)
        image.set_aperture_range(1, 450. * au, 450. * au)
        image.set_wavelength_range(1, 0.9, 1.1)
        image.set_stokes(True)

        m.set_raytracing(True)

        m.set_n_initial_iterations(3)

        m.set_cylindrical_polar_grid_auto(100, 30, 1)

        m.set_monochromatic(True, wavelengths=[1.])

        m.set_mrw(True, gamma=2.)
        # Don't use the PDA here because it's too slow when there are too few photons

        m.set_n_photons(initial=10000, imaging_sources=250000, imaging_dust=500000,
                        raytracing_sources=100000, raytracing_dust=100000)

        m.set_max_interactions(10000000)

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)

    @bit_level
    @pytest.mark.parametrize(('tau'), [1000, 10000, 100000, 1000000])
    def test_pinte_specific_energy(self, tmpdir, tau, generate):

        m = AnalyticalYSOModel()

        # Star:
        #
        # T_star = 4,000K
        # R_star = 2 * R_sun

        m.star.radius = 2. * rsun
        m.star.temperature = 4000.
        m.star.luminosity = 4. * pi * (2. * rsun) ** 2. * sigma * 4000. ** 4.

        # Disk:
        #
        # Sigma(r) = Sigma_0 * (r/r_0)^-1.5
        # h(r) = h_0 * (r/r_0)^1.125
        # h_0 = 10 AU at 100AU
        # r_in = 0.1AU, r_out = 400AU\
        # Sharp vertical edges in cylindrical polar coordinates
        # Disk mass (dust!) = 3e-8, 3e-7, 3e-6, 3e-5

        disk = m.add_flared_disk()
        disk.p = -1.5
        disk.beta = 1.125
        disk.mass = 3.e-8 * msun * tau / 1.e3
        disk.rmin = 0.1 * au
        disk.rmax = 400 * au
        disk.h_0 = 10 * au
        disk.r_0 = 100. * au

        disk.cylindrical_inner_rim = True
        disk.cylindrical_outer_rim = True

        # Dust:
        #
        # Spherical grains, 1 micron, 3.5g/cm^3

        disk.dust = SphericalDust(os.path.join(DATA, 'pinte_dust_lite.hdf5'))

        m.set_n_initial_iterations(3)

        m.set_cylindrical_polar_grid_auto(50, 30, 1)

        m.set_mrw(True, gamma=2.)
        m.set_pda(True)

        m.set_n_photons(initial=50000, imaging=0)

        m.set_max_interactions(10000000)

        m.set_copy_input(False)
        m.write(tmpdir.join(random_id()).strpath, copy=False, absolute_paths=True)
        output_file = tmpdir.join(random_id()).strpath
        m.run(output_file)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".rtout")
            shutil.copy(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".rtout")
            assert_identical_results(output_file, reference_file)
