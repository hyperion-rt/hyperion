Advanced configuration
======================

In :doc:`../setup/setup_conf`, we saw how to set some of the basic parameters
that determine how Hyperion is run, and in this section we discuss some more
advanced options.

.. _monochromatic-rt:

Monochromatic radiative transfer
--------------------------------

By default, when producing images, the radiative transfer is done over all
wavelengths where radiation is emitted. Every emitted photon is propagated,
and may be scattered or absorbed and re-emitted (or both) and we follow the
propagation until the photon escapes the grid. When the photon is binned into
an image or SED (whether making images with the binning or peeling-off
algorithms), the photon is binned into a wavelength grid.

This means that producing images at exact wavelengths can be very inefficient,
because it may be that 99% of the photons end up at a different wavelength and
do not contribute to the images. In order to get around this problem, Hyperion
also implements the concept of monochromatic radiative transfer (see section
2.6.4 of the `Hyperion paper <http://adsabs.harvard.edu/abs/2011A%26A...536A..79R>`_).
In short, the way this algorithm works is the following:

* The temperature is first computed as usual in the initial iterations.

* Photons are then emitted from the source at the wavelengths of interest and
  propagated until they are absorbed, at which point they are terminated. In
  practice, this process can be inefficient in cases where the albedo is low,
  because photons may not often scatter. Therefore, we instead force the photon
  to scatter continuously, but decrease its energy by the albedo at each
  interaction to make sure energy is conserved. Once the energy of the photon
  has decreased below a certain threshold, it is terminated.

* Photons are then emitted from the dust and propagated until they are
  terminated, as for the source photons described above.

* If raytracing is not used, then the photons are peeled off both when the
  photons are emitted from the source or dust, as well as when they scatter.
  Otherwise, they are peeled off only when they scatter, and the source and dust
  components are computed in the raytracing phase.

.. note:: The algorithm that forces scatterings and lowers the energy of the
          photon packets is new in Hyperion v0.9.9. Prior to this, the photons
          were simply propagated until they were absorbed.

To enable monochromatic radiative transfer, before setting the number of
photons, you should call the :meth:`~hyperion.model.Model.set_monochromatic`
method. For example, to compute images at 1, 1.2, and 1.4 microns, you would
need to do::

    m.set_monochromatic(True, wavelengths=[1., 1.2, 1.4])

where the ``wavelength`` arguments takes a list of wavelengths in microns. By
default, the energy threshold for termination is set to 1.e-10 of the original
photon energy, but this can be changed using the ``energy_threshold`` argument::

    m.set_monochromatic(True, wavelengths=[1., 1.2, 1.4],
                        energy_threshold=1.e-20)

When using the monochromatic mode, it is then necessary to set the number of
photons separately for the photons emitted from sources and the photons emitted
from dust::

    m.set_n_photons(..., imaging_sources=1000, imaging_dust=1000, ...)

This should be used instead of the ``imaging`` option. The number of photons
is the number **per wavelength**.

.. _pure-scattering:

Minimum column density to escape each cell
------------------------------------------

In very optically thick models, many photons are emitted or scattered in
regions from which light cannot escape. Hyperion can estimate, for every
cell, the minimum column density (for each dust type) that a photon has to
cross to escape the grid from anywhere in the cell, by casting random rays
through the grid after the initial iterations::

    m.set_compute_escape_column_density(True, rays_per_cell=30)

Every ray provides a sample for all the cells it crosses on the way to the
edge of the grid, so the number of samples per cell is usually much larger
than ``rays_per_cell``. Since the samples are exact column densities along
real paths, the estimate can only be larger than or equal to the true
minimum, and converges to it as the number of rays increases.

The minimum column density is written to the output file, and can be accessed
with :meth:`~hyperion.model.ModelOutput.get_quantities`::

    from hyperion.model import ModelOutput
    m = ModelOutput('model.rtout')
    grid = m.get_quantities()
    column = grid['minimum_escape_column_density'][0].array  # first dust type

Multiplying this by the extinction coefficient (``chi``) of the dust at a
given frequency gives the minimum optical depth to escape each cell at that
frequency.

When the minimum escape column density is computed, photons in cells where
the minimum optical depth to escape is larger than the maximum optical depth
for peeled photons (see
:meth:`~hyperion.model.Model.set_max_peeloff_optical_depth`) are no longer
peeled off, in both normal and monochromatic mode. This avoids integrating
the optical depth for peeled photons whose contribution would be smaller than
``exp(-tau_max)`` (exactly zero in double precision for the default maximum
of 1000).

In monochromatic mode, it is also possible to avoid emitting dust photons
from cells where the minimum optical depth to escape exceeds a given value::

    m.set_compute_escape_column_density(True, rays_per_cell=30)
    m.set_max_peeloff_optical_depth(100.)
    m.set_monochromatic(True, wavelengths=[...], max_emission_escape_optical_depth=200.)

In models where most of the dust emission at a given wavelength comes from
regions that are very optically thick at that wavelength (for example the
inner midplane of a massive disk at near- and mid-infrared wavelengths), this
can greatly improve the signal-to-noise of the scattered dust emission for a
given computation time, since photons are then only emitted from the
regions they can escape from. This removes the contribution from the other
cells, so the value should be chosen to be sufficiently larger than the
maximum optical depth for peeled photons, since photons emitted where the
optical depth to escape is larger have to scatter many times through very
optically thick regions before being able to contribute. This requires the
minimum escape column density to be computed, and since the minimum optical
depth is computed to the edge of the grid, this cannot be used with observers
inside the grid or with peeled images that ignore the optical depth (an
error is raised when the model is written in all these cases).

Scattered-light images
----------------------

In some cases, one might want to compute scattered light images at wavelengths
where there is no dust emission. In this case, there is no need to compute the
specific energy of the dust, and there is also no need in re-emitting photons
when computing images/SEDs. Therefore, one can set::

    m.set_n_initial_iterations(0)
    m.set_kill_on_absorb(True)
    m.set_raytracing(True)

which turns off the specific energy calculation, kills photons as soon as they
are first absorbed, and enables raytracing for the source emission. For the
photon numbers, one can set ``raytracing_dust=0`` to zero, since this is not
needed (there is no dust emission).

.. note:: This cannot be used for *all* scattered light images. For example,
          in a protostar, a K-band image may have a non-negligeable amount of
          scattered light flux originating from the inner rim of the disk.
          This technique can only be used when there is no dust emission.

This can be combined with the `Monochromatic radiative transfer`_ option
described above to avoid wasting photons at wavelengths where they are not
needed. When treating only scattering, you will then want to set the following
options::

    m.set_n_photons(imaging_sources=1000, imaging_dust=0,
                    raytracing_sources=1000, raytracing_dust=0)

where the values should be adjusted to your model, but the important point is
that ``initial`` is not needed, and ``imaging_dust`` and ``raytracing_dust``
can be set to 0.

For a full example of a model computing scattered light images, see
:doc:`../tutorials/howto_pure_scattering`.

Miscellaneous Settings
----------------------

Set the maximum number of photon interactions::

    m.set_max_interactions(100000)

When computing the contribution of photons to peeled images and SEDs, the
optical depth along the line of sight is only integrated up to an optical depth
of 1000, beyond which the contribution of the photon is exactly zero in double
precision. In very optically thick models, the imaging can be made much faster
by lowering this value with
:meth:`~hyperion.model.Model.set_max_peeloff_optical_depth`::

    m.set_max_peeloff_optical_depth(100.)

Contributions attenuated by more than exp(-100) are then neglected, which is
only noticeable in parts of images that are fainter than the brightest parts
by many orders of magnitude.

Set the number of output bytes per floating point value for the physical
arrays (4 = 32-bit, 8 = 64-bit)::

    m.set_output_bytes(4)

To set the minimum temperature for dust::

    m.set_minimum_temperature(10.)
    m.set_minimum_temperature([10., 5., 20.])

If a scalar value is specified, the same value is used for all dust types. If
a list is specified, the list should have as many items as dust types, and
each item corresponds to the minimum temperature for each dust type.

Similarly, to set the minimum specific energy::

    m.set_minimum_specific_energy(1.e-4)
    m.set_minimum_specific_energy([1.e-4, 1.e-5, 2.e-5])

By default, photon positions and cells are double-checked every 1 in 1000 cell
crossings. This can be changed
with :meth:`~hyperion.model.Model.set_propagation_check_frequency`::

    m.set_propagation_check_frequency(0.01)

Note that values higher than 0.001 (the default) will cause the code to slow
down. Setting the frequency to zero disables the checks entirely, which makes
the code slightly faster (by a few percent) since no random number is then
drawn at each cell crossing::

    m.set_propagation_check_frequency(0)
