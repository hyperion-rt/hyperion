Important Notes
===============

Gridding
--------

As a user, you are responsible for ensuring that the grid on which the
radiative transfer is carried out is adequate to resolve density and
temperature gradients. In the case of the
:class:`~hyperion.model.AnalyticalYSOModel` class, Hyperion tries to optimize
the gridding of the inner disk to ensure that it is properly resolved, but you
still need to ensure that the *number* of cells you specificy is adequate.


Number of photons
-----------------

Similarly to the `Gridding`_, you are responsible for choosing an adequate
number of photons to carry out the radiative transfer. You can read the
:doc:`../setup/photon_numbers` page for advice, and this page will be expanded
in future, but the best way to ensure that you have enough photons is to check
yourself that the results (temperature, images, SEDs) have converged if you
increase the number of photons.

Dust and/or Gas?
----------------

While Hyperion is a dust continuum radiative transfer code, you should take
care when specifying densities, accretion rates, etc. as to whether to include
the contribution from gas. The guidelines are as follows:

* The core Fortran code does not make any assumptions regarding gas - it simply
  computes optical depths from the densities and opacities provided. Therefore,
  dust opacities and densities have to be consistent as to whether they are per
  unit dust mass, or per unit dust+gas mass, and what gas-to-dust ratio they
  assume.

  For example, if the dust opacities are provided per unit dust
  mass, then the densities specified in the model should be dust densities. If
  the dust opacities are provided per unit dust+gas mass, then the densities
  specified in the model should be dust+gas densities. For the dust models
  provided in :doc:`../dust/dust`, each dust model explicitly states whether it
  includes gas in the opacities or not. When setting up your own dust models,
  you should be aware of whether they are given per unit dust or dust+gas mass
  and what gas-to-dust ratio was assumed.

* For the :class:`~hyperion.densities.UlrichEnvelope` density structure, the
  infall rate provided is directly related to the density, so that if the dust
  opacities are per unit dust mass, the infall rate should be the infall rate
  of dust.

* For the :class:`~hyperion.densities.AlphaDisk` density structure, the
  accretion rate provided does not relate to the density in the model, but it
  used to add a source of luminosity. Therefore, it should include the total
  mass of **dust+gas** regardless of how the opacities are expressed. However,
  the mass or density of the disk should be given depending on the units of the
  opacity, as described above.

