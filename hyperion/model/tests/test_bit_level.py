# The purpose of these tests is to check that the results produced by
# Hyperion agree at the bit-level. The data for these tests will need to be
# changed if algorithms are updated, and bugs are fixed, but otherwise there
# is no reason we should expect a change in results from one commit to the
# next. Since these files take place, we should minimize the number of tests
# to run.

import os
import itertools

import h5py
import pytest
import numpy as np

import cPickle as pickle

from .test_helpers import random_filename
from .. import Model, AnalyticalYSOModel
from ...util.constants import pc, lsun, c, au, msun, pi, sigma, rsun
from ...grid import CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid, AMRGrid, OctreeGrid
from ...dust import IsotropicDust

GRID_TYPES = ['car', 'cyl', 'sph', 'amr', 'oct']

DATA = os.path.join(os.path.dirname(__file__), 'data')

generate_reference = pytest.mark.generate_reference


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
    level = self.grid['amr'].add_level()
    grid = level.add_grid()
    grid.xmin, grid.xmax = -u, u
    grid.ymin, grid.ymax = -u, u
    grid.zmin, grid.zmax = -u, u
    grid.nx, grid.ny, grid.nz = 8, 6, 4
    grid.quantities['density'] = np.random.random((4, 6, 8)) * d
    grid.quantities['density_2'] = np.random.random((4, 6, 8)) * d
    grid.quantities['density_3'] = np.random.random((4, 6, 8)) * d

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
        if arg not in ['self', 'generate']:
            name += ["{0}={1}".format(arg, values[arg])]
    name = '.'.join(name)
    return name


def assert_output_matches(filename, reference):

    differences = []

    # Get items from file to compare
    groups, datasets, attributes = make_item_list(filename)

    # Read in reference from pickle
    ref = open(reference, 'rb')
    groups_ref = pickle.load(ref)
    datasets_ref = pickle.load(ref)
    attributes_ref = pickle.load(ref)

    # Order groups since order is not guaranteed
    groups.sort()
    groups_ref.sort()

    # Check that group lists are the same
    if groups != groups_ref:
        differences.append("Group lists do not match: found {0} but expected {1}".format(str(groups), str(groups_ref)))

    # Make ordered lists of the datasets to compare
    dataset_list = sorted(datasets.keys())
    dataset_ref_list = sorted(datasets_ref.keys())

    # Check whether the dataset lists are different
    if dataset_list != dataset_ref_list:
        differences.append("Dataset lists do not match: found {0} but expected {1}".format(str(dataset_list), str(dataset_ref_list)))
    else: # Check that hashes match
        for d in datasets:
            if datasets[d] != datasets_ref[d]:
                differences.append("Dataset hashes do not match: found {0}={1} but expected {0}={2}".format(d, datasets[d], datasets_ref[d]))

    # Make ordered lists of the attributes to compare
    attribute_list = sorted(attributes.keys())
    attribute_ref_list = sorted(attributes_ref.keys())

    # Check whether the attribute lists are different
    if attribute_list != attribute_ref_list:
        differences.append("Attribute lists do not match: found {0} but expected {1}".format(str(attribute_list), str(attribute_ref_list)))
    else: # Check that hashes match
        for a in attributes:
            if attributes[a] != attributes_ref[a]:
                differences.append("Attribute values do not match: found {0}={1} but expected {0}={2}".format(a, attributes[a], attributes_ref[a]))

    for item in differences:
        print(item)

    assert len(differences) == 0


def write_item_list(filename, filename_out):

    groups, datasets, attributes = make_item_list(filename)

    f = open(filename_out, 'wb')
    pickle.dump(groups, f, 2)
    pickle.dump(datasets, f, 2)
    pickle.dump(attributes, f, 2)
    f.close()


def type_cast(a):
    try:
        float(a)
        if float(a) == int(a):
            return int(a)
        else:
            return float(a)
    except:
        try:  # in case it is a bytes object
            return a.decode()
        except:
            return str(a)


def make_item_list(filename):

    from hashlib import md5

    # List of attributes to exclude from checking (time-dependent)
    EXCLUDE_ATTR = ['date_started', 'date_ended', 'cpu_time', 'python_version', 'fortran_version']

    groups = []
    datasets = {}
    attributes = {}

    # Open file
    f = h5py.File(filename, 'r')

    # List datasets and groups in file
    def func(name, obj):
        if isinstance(obj, h5py.Dataset):
            a = np.array(obj)
            datasets[name] = md5(a).hexdigest()
        elif isinstance(obj, h5py.Group):
            groups.append(name)
    f.visititems(func)

    # Loop over all groups and datasets to check attributes
    for item in ['/'] + datasets.keys() + groups:

        # Find all attributes
        attr = f[item].attrs.keys()

        for a in attr:
            if a not in EXCLUDE_ATTR and a.lower() == a:
                attributes[a] = type_cast(f[item].attrs[a])

    return groups, datasets, attributes


class TestEnergy(object):

    def setup_class(self):
        setup_all_grid_types(self, pc, 1.e-20)
        self.dust_file = os.path.join(DATA, 'kmh_lite.hdf5')

    @generate_reference
    @pytest.mark.parametrize(('grid_type', 'sample_sources_evenly', 'multiple_densities'), list(itertools.product(GRID_TYPES, [False, True], [False, True])))
    def test_specific_energy(self, grid_type, sample_sources_evenly, multiple_densities, generate=False):

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

        # m.set_copy_input(False)

        m.set_sample_sources_evenly(sample_sources_evenly)

        m.conf.output.output_specific_energy = 'all'

        m.write(random_filename())
        output_file = random_filename()
        m.run(output_file, overwrite=True)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".pickle")
            write_item_list(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".pickle")
            assert_output_matches(output_file, reference_file)

    @generate_reference
    @pytest.mark.parametrize(('grid_type', 'raytracing', 'sample_sources_evenly'), list(itertools.product(GRID_TYPES, [False, True], [False, True])))
    def test_peeloff(self, grid_type, raytracing, sample_sources_evenly, generate=False):

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

        # m.set_copy_input(False)

        m.set_sample_sources_evenly(sample_sources_evenly)

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(5, 0.05, 200.)
        i_p.set_viewing_angles([33.4, 110.], [65.4, 103.2])
        i_p.set_image_size(4, 5)
        i_p.set_image_limits(-0.8 * pc, 0.8 * pc, -pc, pc)
        i_p.set_aperture_range(5, 0.1 * pc, pc)

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(4, 0.05, 200.)
        i_p.set_viewing_angles([22.1], [203.2])
        i_p.set_image_size(6, 6)
        i_p.set_image_limits(-pc, pc, -pc, pc)
        i_p.set_aperture_range(2, 0.5 * pc, pc)
        i_p.set_track_origin('basic')

        i_p = m.add_peeled_images()
        i_p.set_wavelength_range(4, 0.05, 200.)
        i_p.set_viewing_angles([22.1], [203.2])
        i_p.set_image_size(6, 6)
        i_p.set_image_limits(-pc, pc, -pc, pc)
        i_p.set_aperture_range(2, 0.5 * pc, pc)
        i_p.set_track_origin('detailed')

        m.write(random_filename())
        output_file = random_filename()
        m.run(output_file, overwrite=True)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".pickle")
            write_item_list(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".pickle")
            assert_output_matches(output_file, reference_file)


class TestPascucci(object):
    """
    Run a very low signal-to-noise version of the Pascucci benchmark models.
    The point is not to have high signal-to-noise but to look for bit-level
    changes.
    """

    def setup_class(self):

        optSi = """
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
            """

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
        self.dust = IsotropicDust(nu[::-1], albedo[::-1], chi[::-1])
        self.dust.optical_properties._extrapolate(1.e-3, 1.e5)
        self.dust.set_lte_emissivities(n_temp=100, temp_min=0.1, temp_max=1600.)

    @generate_reference
    @pytest.mark.parametrize(('tau'), [0.1, 1, 10, 100])
    def test_pascucci(self, tau, generate=False):

        m = AnalyticalYSOModel()

        # Star:
        #
        # T_star = 5,800K
        # R_star = R_sun
        # M_star = M_sun

        m.star.radius = 1. # = 1 cm (i.e. point source)
        m.star.temperature = 5800.
        m.star.luminosity = 4. * pi * rsun**2 * sigma * 5800.**4

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

        disk.dust = self.dust

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
        image.set_image_limits(-1500.*au, 1500.*au, -1500.*au, 1500.*au)
        image.set_aperture_range(1, 1500.*au, 1500.*au)
        image.set_wavelength_range(61, 1, 61)

        m.set_raytracing(True)
        m.set_n_initial_iterations(5)
        # Don't test for convergence, since this is a low signal-to-noise models

        # Use a lower-resolution grid
        m.set_spherical_polar_grid_auto(100, 30, 1, rmax=1300.*au)

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

        m.write(random_filename())
        output_file = random_filename()
        m.run(output_file, overwrite=True)

        if generate:
            reference_file = os.path.join(generate, function_name() + ".pickle")
            write_item_list(output_file, reference_file)
            pytest.skip("Skipping test, since generating data")
        else:
            reference_file = os.path.join(DATA, function_name() + ".pickle")
            assert_output_matches(output_file, reference_file)
