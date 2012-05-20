# The purpose of these tests is to ensure that photons are not getting
# killed during propagation due to special alignments or source positions.
# We don't use parametrized tests here, because each grid has different
# special alignemnts.

import h5py
import pytest
import numpy as np
from .. import Model
from .test_helpers import random_filename, get_test_dust


def any_photons_killed(handle):
    return handle.attrs['killed_photons_geo_initial'] != 0 or \
           handle.attrs['killed_photons_int_initial'] != 0 or \
           handle.attrs['killed_photons_geo_final'] != 0 or \
           handle.attrs['killed_photons_int_final'] != 0 or \
           handle.attrs['killed_photons_geo_raytracing'] != 0 or \
           handle.attrs['killed_photons_int_raytracing'] != 0

# SPHERICAL POLAR GRID


class TestSphericalBase(object):

    r = np.linspace(0., 10., 15)
    t = np.linspace(0., np.pi, 17)
    p = np.linspace(0., 2. * np.pi, 15)

    def setup_method(self, method):

        dust = get_test_dust()

        self.m = Model()
        self.m.set_spherical_polar_grid(self.r, self.t, self.p)
        self.m.add_density_grid(np.ones(self.m.grid.shape) * 1.e-20, dust)

    def test_ptsource_origin(self):
        '''A point source exactly at the origin'''

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(random_filename())
        file_out = random_filename()
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_vertices(self):
        '''Point sources exactly on the vertices'''

        for ir in range(len(self.r) - 1):
            for it in range(len(self.t)):
                for ip in range(len(self.p)):
                    s = self.m.add_point_source()
                    x = self.r[ir] * np.cos(self.p[ip]) * np.sin(self.t[it])
                    y = self.r[ir] * np.sin(self.p[ip]) * np.sin(self.t[it])
                    z = self.r[ir] * np.cos(self.t[it])
                    # Clip to w=0 (necessary due to numerical precision)
                    x = 0 if abs(x) < 1.e-10 else x
                    y = 0 if abs(y) < 1.e-10 else y
                    s.position = (x, y, z)
                    s.luminosity = 1
                    s.temperature = 5000.

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(random_filename())
        file_out = random_filename()
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    # The following test is known to fail - that is, if sources are very
    # close to, but not on w=0, then photons are killed.
    @pytest.mark.xfail()
    def test_ptsource_vertices_noclip(self):
        '''Point sources on vertices, but without clipping'''

        for ir in range(len(self.r) - 1):
            for it in range(len(self.t)):
                for ip in range(len(self.p)):
                    s = self.m.add_point_source()
                    x = self.r[ir] * np.cos(self.p[ip]) * np.sin(self.t[it])
                    y = self.r[ir] * np.sin(self.p[ip]) * np.sin(self.t[it])
                    z = self.r[ir] * np.cos(self.t[it])
                    s.position = (x, y, z)
                    s.luminosity = 1
                    s.temperature = 5000.

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(random_filename())
        file_out = random_filename()
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff(self):
        '''Peeloff angles aligned with walls for source at origin'''

        theta = np.linspace(0., 180., 37)
        phi = np.linspace(0., 360., 73)

        T, P = np.meshgrid(theta, phi)

        theta = T.ravel()
        phi = P.ravel()

        i = self.m.add_peeled_images()
        i.set_viewing_angles(theta, phi)
        i.set_image_size(1, 1)
        i.set_image_limits(-self.r[-1], self.r[-1], -self.r[-1], self.r[-1])
        i.set_wavelength_range(1, 0.1, 10.)

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(0)
        self.m.set_n_photons(imaging=100)

        self.m.write(random_filename())
        file_out = random_filename()
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)
