
import os
import tempfile

import numpy as np
import pytest

from .. import Model
from ...util.functions import random_id


@pytest.mark.requires_hyperion_binaries
def test_spot_uses_its_own_spectrum():
    # Regression test: photons emitted from a spot must use the spot's own
    # spectrum, not the parent sphere's. We give the sphere and the spot
    # disjoint emission bands (no dust, so photons escape directly), so any
    # flux in the spot's band can only come from the spot having used its own
    # spectrum. Before the fix, source_emit overwrote the spot-sampled
    # frequency with the sphere spectrum, so the spot band got zero flux.

    m = Model()
    m.set_cartesian_grid([-1e12, 1e12], [-1e12, 1e12], [-1e12, 1e12])

    nu = np.logspace(np.log10(3e12), np.log10(1e15), 300)
    fnu_sphere = np.where((nu > 1e13) & (nu < 2e13), 1., 0.)   # ~15-30 micron
    fnu_spot = np.where((nu > 3e14) & (nu < 6e14), 1., 0.)     # ~0.5-1 micron (disjoint)

    s = m.add_spherical_source()
    s.radius = 1e11
    s.position = (0., 0., 0.)
    s.luminosity = 1.
    s.spectrum = (nu, fnu_sphere)

    spot = s.add_spot()
    spot.longitude = 0.
    spot.latitude = 0.
    spot.radius = s.radius
    spot.luminosity = 1.
    spot.spectrum = (nu, fnu_spot)

    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_viewing_angles([45.], [45.])
    sed.set_wavelength_range(60, 0.1, 100.)
    sed.set_aperture_radii(1, 1e12, 1e12)

    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=20000)

    tmpdir = tempfile.mkdtemp()
    m.write(os.path.join(tmpdir, random_id()))
    mo = m.run()

    wav, nufnu = mo.get_sed()
    wav = np.array(wav)
    nufnu = np.squeeze(np.array(nufnu))

    sphere_band = (wav > 10.) & (wav < 40.)    # ~15-30 micron
    spot_band = (wav > 0.4) & (wav < 1.2)      # ~0.5-1 micron

    sphere_flux = np.nansum(nufnu[sphere_band])
    spot_flux = np.nansum(nufnu[spot_band])

    # the model actually ran and the sphere emitted in its band
    assert sphere_flux > 0.
    # the spot emitted in its OWN band (exactly zero before the fix)
    assert spot_flux > 0.
    assert spot_flux > 0.2 * sphere_flux
    # and essentially nothing leaked between the two disjoint bands
    assert spot_flux + sphere_flux > 0.9 * np.nansum(nufnu)
