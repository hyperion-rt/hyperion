Advanced configuration
======================

In :doc:`../setup/setup_conf`, we saw how to set some of the basic parameters that determine how Hyperion is run, and in this section we discuss some more advanced options.

Monochromatic radiative transfer
--------------------------------

.. note:: This section is being written

Scattered-light images
----------------------

In some cases, one might want to compute scattered light images at wavelengths where there is no dust emission. In this case, there is no need to compute the specific energy of the dust, and there is also no need in re-emitting photons when computing images/SEDs. Therefore, one can set::

    m.set_n_initial_iterations(0)
    m.set_kill_on_absorb(True)
    m.set_raytracing(True)

which turns off the specific energy calculation, kills photons as soon as they are first absorbed, and enables raytracing for the source emission. For the photon numbers, one can set ``raytracing_dust=0`` to zero, since this is not needed (there is no dust emission).

.. note:: This cannot be used for *all* scattered light images. For example,
          in a protostar, a K-band image may have a non-negligeable amount of
          scattered light flux originating from the inner rim of the disk.
          This technique can only be used when there is no dust emission.

Miscellaneous Settings
----------------------

Set the maximum number of photon interactions::

    m.set_max_interactions(100000)

Set the number of output bytes per floating point value for the physical
arrays (4 = 32-bit, 8 = 64-bit)::

    m.set_output_bytes(4)

For the :doc:`model` class, one can set a minimum temperature to which
temperatures below this will be reset::

    m.add_density_grid(density, dust, minimum_temperature=100.)

and in terms of specific energy::

    m.add_density_grid(density, dust, minimum_specific_energy=100.)

For the :doc:`analytical_yso_model`, this can be done with::

    m.set_minimum_temperature(100.)