# The purpose of these tests is to ensure that photons are not getting
# killed during propagation due to special alignments or source positions.
# We don't use parametrized tests here, because each grid has different
# special alignemnts.

import h5py
from astropy.tests.helper import pytest
import numpy as np
from .. import Model
from .test_helpers import random_id, get_test_dust


def any_photons_killed(handle):
    for group in handle:
        if 'killed_photons_geo' in handle[group].attrs:
            if handle[group].attrs['killed_photons_geo'] != 0 or handle[group].attrs['killed_photons_int'] != 0:
                return True
    return handle.attrs['killed_photons_geo_final'] != 0 or \
           handle.attrs['killed_photons_int_final'] != 0 or \
           handle.attrs['killed_photons_geo_raytracing'] != 0 or \
           handle.attrs['killed_photons_int_raytracing'] != 0


class TestCartesianBase(object):

    x = np.linspace(-10., 10., 15)
    y = np.linspace(-10., 10., 15)
    z = np.linspace(-10., 10., 15)

    def setup_method(self, method):

        dust = get_test_dust()

        self.m = Model()
        self.m.set_cartesian_grid(self.x, self.y, self.z)
        self.m.add_density_grid(np.ones(self.m.grid.shape) * 1.e-40, dust)

    def test_ptsource_origin(self, tmpdir):
        '''A point source exactly at the origin'''

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_vertices(self, tmpdir):
        '''Point sources exactly on the vertices'''

        for ix in range(1, len(self.x) - 1):
            for iy in range(1, len(self.y) - 1):
                for iz in range(1, len(self.z) - 1):
                    s = self.m.add_point_source()
                    x = self.x[ix]
                    y = self.y[iy]
                    z = self.z[iz]
                    s.position = (x, y, z)
                    s.luminosity = 1
                    s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    @pytest.mark.xfail()
    def test_ptsource_edge(self, tmpdir):
        '''Point sources exactly on the edge of the grid'''

        for ix in range(len(self.x)):
            for iy in range(len(self.y)):
                for iz in range(len(self.z)):
                    if ix == 0 or ix == len(self.x) - 1 or \
                        iy == 0 or iy == len(self.y) - 1 or \
                            iz == 0 or iz == len(self.z) - 1:
                        s = self.m.add_point_source()
                        x = self.x[ix]
                        y = self.y[iy]
                        z = self.z[iz]
                        s.position = (x, y, z)
                        s.luminosity = 1
                        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff_aligned(self, tmpdir):
        '''Peeloff for source at origin'''

        theta = np.linspace(0., 180., 37)
        phi = np.linspace(0., 360., 73)

        T, P = np.meshgrid(theta, phi)

        theta = T.ravel()
        phi = P.ravel()

        i = self.m.add_peeled_images()
        i.set_viewing_angles(theta, phi)
        i.set_image_size(1, 1)
        i.set_image_limits(-self.x[-1], self.x[-1], -self.x[-1], self.x[-1])
        i.set_wavelength_range(1, 0.1, 10.)

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(0)
        self.m.set_n_photons(imaging=100)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)


class TestCartesianLarge(TestCartesianBase):

    x = np.linspace(-10.e20, 10.e20, 15)
    y = np.linspace(-10.e20, 10.e20, 15)
    z = np.linspace(-10.e20, 10.e20, 15)


class TestCartesianSmall(TestCartesianBase):

    x = np.linspace(-10.e-20, 10.e-20, 15)
    y = np.linspace(-10.e-20, 10.e-20, 15)
    z = np.linspace(-10.e-20, 10.e-20, 15)


class TestSphericalBase(object):

    r = np.linspace(0., 10., 15)
    t = np.linspace(0., np.pi, 17)
    p = np.linspace(0., 2. * np.pi, 15)

    def setup_method(self, method):

        dust = get_test_dust()

        self.m = Model()
        self.m.set_spherical_polar_grid(self.r, self.t, self.p)
        self.m.add_density_grid(np.ones(self.m.grid.shape) * 1.e-40, dust)

    def test_ptsource_origin(self, tmpdir):
        '''A point source exactly at the origin'''

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_vertices(self, tmpdir):
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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    @pytest.mark.xfail()
    def test_ptsource_edge(self, tmpdir):
        '''Point sources exactly on the vertices'''

        ir = -1
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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    # The following test is known to fail - that is, if sources are very
    # close to, but not on w=0, then photons are killed.
    @pytest.mark.xfail()
    def test_ptsource_vertices_noclip(self, tmpdir):
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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff(self, tmpdir):
        '''Peeloff for source at origin'''

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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff_aligned(self, tmpdir):
        '''Peeloff for source at origin'''

        theta = self.t
        phi = self.p

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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)


class TestSphericalLarge(object):

    r = np.linspace(0., 10.e20, 15)
    t = np.linspace(0., np.pi, 17)
    p = np.linspace(0., 2. * np.pi, 15)


class TestSphericalSmall(object):

    r = np.linspace(0., 10.e-20, 15)
    t = np.linspace(0., np.pi, 17)
    p = np.linspace(0., 2. * np.pi, 15)


class TestCylindricalBase(object):

    w = np.linspace(0., 10., 15)
    z = np.linspace(-5., 5., 17)
    p = np.linspace(0., 2. * np.pi, 15)

    def setup_method(self, method):

        dust = get_test_dust()

        self.m = Model()
        self.m.set_cylindrical_polar_grid(self.w, self.z, self.p)
        self.m.add_density_grid(np.ones(self.m.grid.shape) * 1.e-40, dust)

    def test_ptsource_origin(self, tmpdir):
        '''A point source exactly at the origin'''

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_vertices(self, tmpdir):
        '''Point sources exactly on the vertices'''

        for iw in range(len(self.w) - 1):
            for iz in range(len(self.z)):
                for ip in range(len(self.p)):
                    s = self.m.add_point_source()
                    x = self.w[iw] * np.cos(self.p[ip])
                    y = self.w[iw] * np.sin(self.p[ip])
                    z = self.z[iz]
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

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    @pytest.mark.xfail()
    def test_ptsource_edge(self, tmpdir):
        '''Point sources exactly on the vertices'''

        iw = -1
        for iz in [self.z[0], self.z[-1]]:
            for ip in range(len(self.p)):
                s = self.m.add_point_source()
                x = self.w[iw] * np.cos(self.p[ip])
                y = self.w[iw] * np.sin(self.p[ip])
                z = self.z[iz]
                # Clip to w=0 (necessary due to numerical precision)
                x = 0 if abs(x) < 1.e-10 else x
                y = 0 if abs(y) < 1.e-10 else y
                s.position = (x, y, z)
                s.luminosity = 1
                s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_vertices_noclip(self, tmpdir):
        '''Point sources on vertices, but without clipping'''

        for iw in range(len(self.w) - 1):
            for iz in range(len(self.z)):
                for ip in range(len(self.p)):
                    s = self.m.add_point_source()
                    x = self.w[iw] * np.cos(self.p[ip])
                    y = self.w[iw] * np.sin(self.p[ip])
                    z = self.z[iz]
                    s.position = (x, y, z)
                    s.luminosity = 1
                    s.temperature = 5000.

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(1)
        self.m.set_n_photons(initial=100000, imaging=0)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff(self, tmpdir):
        '''Peeloff for source at origin'''

        theta = np.linspace(0., 180., 37)
        phi = np.linspace(0., 360., 73)

        T, P = np.meshgrid(theta, phi)

        theta = T.ravel()
        phi = P.ravel()

        i = self.m.add_peeled_images()
        i.set_viewing_angles(theta, phi)
        i.set_image_size(1, 1)
        i.set_image_limits(-self.w[-1], self.w[-1], -self.w[-1], self.w[-1])
        i.set_wavelength_range(1, 0.1, 10.)

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(0)
        self.m.set_n_photons(imaging=100)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)

    def test_ptsource_origin_peeloff_aligned(self, tmpdir):
        '''Peeloff for source at origin'''

        theta = np.linspace(0., 180., 37)
        phi = self.p

        T, P = np.meshgrid(theta, phi)

        theta = T.ravel()
        phi = P.ravel()

        i = self.m.add_peeled_images()
        i.set_viewing_angles(theta, phi)
        i.set_image_size(1, 1)
        i.set_image_limits(-self.w[-1], self.w[-1], -self.w[-1], self.w[-1])
        i.set_wavelength_range(1, 0.1, 10.)

        s = self.m.add_point_source()
        s.position = (0., 0., 0.)
        s.luminosity = 1
        s.temperature = 5000.

        self.m.set_n_initial_iterations(0)
        self.m.set_n_photons(imaging=100)

        self.m.write(tmpdir.join(random_id()).strpath)
        file_out = tmpdir.join(random_id()).strpath
        self.m.run(file_out)
        f = h5py.File(file_out)
        assert not any_photons_killed(f)


class TestCylindricalLarge(object):

    w = np.linspace(0., 10.e20, 15)
    z = np.linspace(-5.e20, 5.e20, 17)
    p = np.linspace(0., 2. * np.pi, 15)


class TestCylindricalSmall(object):

    w = np.linspace(0., 10.e-20, 15)
    z = np.linspace(-5.e-20, 5.e-20, 17)
    p = np.linspace(0., 2. * np.pi, 15)
