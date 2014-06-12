Moving sources/dust
===================

.. note:: This feature is experimental, and should be used with caution!

About
-----

Hyperion has basic support for computing the Doppler shift in the emission due
to moving sources and moving dust. Of course, in most cases this makes no
difference if both the emission spectrum from the source and the
opacity/emissivities for the dust are relatively featureless (at least on the
frequency scale for which Doppler shifts apply). This feature was implemented
for the case where one includes spectral lines in the spectrum from the star,
and wants to correctly track the strength and broadening of these lines due to
the motion of the stars themselves, and the scattering on the dust if it is
moving.

In detail, the algorithm:

- Emits photons in the frame of reference of the star, and Doppler shifts them
  into the frame of reference of the radiative transfer along the propagation
  direction.

- When scattering, transforms the frequency into the frame of reference of the
  dust, then back into the frame of reference of the radiative transfer, along
  the new propagation direction.

- When emitting from dust, emits photons in the frame of reference of the dust,
  and transforms the frequency into the frame of reference of the radiative
  transfer, along the new propagation direction.

Currently, the algorithm makes the following approximations:

- The opacity of the dust used to compute the optical depth through cells is
  computed from the frequency in the frame of reference of the radiative
  transfer, not the dust, because using the rest frequency would have a
  significant impact computationally. However, the opacity of dust usually
  never changes abruptly enough for this to matter.

- The Doppler shift calculation does take into account the relativistic term,
  but we do not calculate the full Lorentz boost, which changes not only the
  frequency but the direction of propagation of the photons.

- Raytracing cannot be used in the relativistic mode.

These assumptions/approximations are not intended to be permanent, and in
future we should be able to address them all (if there is a need for it).

Setting velocities for sources
------------------------------

All sources support setting the ``velocity`` parameter, which can be set to an
iterable of three values ``(vx, vy, vz)`` similarly to the ``position``
parameter::

    s = m.add_spherical_source()
    s.radius = rsun
    s.luminosity = lsun
    s.spectrum = (nu, fnu)
    s.velocity = [2e8, 0., 0.]

The velocity should be given in c.g.s (``cm/s``).

Setting velocities for dust
---------------------------

To set the velocities of the dust, use the ``velocity`` keyword argument in
``add_density_grid``::

    m.add_density_grid(rho, 'gray_dust.hdf5', velocity=(vx, vy, vz))
    
where ``vx``, ``vy``, and ``vz`` should be objects similar to the density and
specific energy (if specified), containing the velocity values in c.g.s
(``cm/s``). For example, for the 3-d regular grids (cartesian and polar), these
can be 3-d Numpy arrays.
