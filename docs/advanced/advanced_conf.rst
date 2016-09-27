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
In short, the way this algorithm works is that since the temperature has
already been computed by the time the images are being computed, it is
possible to consider only the propagation, scattering, and absorption of
photons at the specific wavelengths/frequencies of interest.

To enable monochromatic radiative transfer, before setting the number of
photons, you should call the :meth:`~hyperion.model.Model.set_monochromatic`
method. For example, to compute images at 1, 1.2, and 1.4 microns, you would need to do::

    m.set_monochromatic(True, wavelengths=[1., 1.2, 1.4])

where the ``wavelength`` arguments takes a list of wavelengths in microns. When
using the monochromatic mode, it is then necessary to set the number of photons
separately for the photons emitted from sources and the photons emitted from
dust::

    m.set_n_photons(..., imaging_sources=1000, imaging_dust=1000, ...)

This should be used instead of the ``imaging`` option. The number of photons
is the number **per wavelength**.

.. _pure-scattering:

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
down.
