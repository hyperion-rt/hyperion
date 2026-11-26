=================
Dust HDF5 Format
=================

Overview
========

Hyperion requires information on the dust properties, such as the albedo,
opacity, mean opacities, and emissivities. These need to be packaged in an
HDF5 format that is described below. In most cases you do not need to create
these files from scratch, and you can instead use the Hyperion Python library
to produce these (see :doc:`../setup/setup_dust`). If however you want to
write the files directly without the Python library, this section is for you!

.. _specification:

Dust file HDF5 format specification
===================================

An HDF5 dust file should contain 5 datasets. The root of the file should contain the following attributes:

* ``emissvar``: whether the emissivity is specified as a function of specific
  energy absorbed in each cell (``E``) or another quantity (but this is not
  supported at this time).

* ``version``: this should be set to ``2`` - the version described in this
  section. Files with version ``1`` did not include the reciprocal Planck
  mean opacities and contained incorrect Rosseland mean opacities, and
  cannot be used for models that use the partial diffusion approximation.

* ``type``: the dust file type. At the moment, the only kind of dust file
  supported is one giving the four unique elements of the scattering matrix
  of dust as a function of scattering angle (``1``). In future, other types
  of dust, such as aligned grains, which require the full 16 elements, will
  be implemented.

* ``lte``: whether the dust emissivities assume local thermodynamic
  equilibrium (LTE).

* ``python_version``: the version of the Python Hyperion library used to
  generate the file. Set this to '0.8.7' if you are writing files yourself
  rather than using the Hyperion library.

The datasets present should be the following:

``optical_properties``
----------------------

This dataset should consist of a table with the basic optical properties of
the dust as a function of frequency, in a binary table. The columns should be:

* ``nu``: The frequency, in Hz.

* ``albedo``: The albedo of the dust.

* ``chi``: The opacity of the dust to extinction, in cm^2/g.

* ``P1``, ``P2``, ``P3``, and ``P4``: The four elements of the scattering
  matrix. These columns should be vector columns, with each table cell
  containing the elements for as many angles as specified in
  ``scattering_angles``.

``scattering_angles``
---------------------

This dataset should consist of a single-column table. The column should be
``mu``, and should give the values of the cosine of the scattering angle for
which the matrix elements are tabulated in the ``Optical properties`` dataset.

``mean_opacities``
------------------

This dataset should consist of a table with pre-computed mean opacities for
the dust. The columns should be:

* ``temperature``: The temperature for which the mean opacities are given,
  in K. This column is not used by the Fortran code, but is needed to read
  the file with the Python library.

* ``specific_energy``: The specific energy corresponding to each temperature,
  in erg/s/g, which is given by :math:`4\sigma T^4\kappa_{\rm P}` where
  :math:`\kappa_{\rm P}` is the Planck mean opacity to absorption.

* ``chi_planck``: The Planck mean opacity to extinction, in cm^2/g

* ``kappa_planck``: The Planck mean opacity to absorption, in cm^2/g

* ``chi_inv_planck``: The reciprocal Planck mean opacity to extinction, in
  cm^2/g

* ``kappa_inv_planck``: The reciprocal Planck mean opacity to absorption, in
  cm^2/g

* ``chi_rosseland``: The Rosseland mean opacity to extinction, in cm^2/g

* ``kappa_rosseland``: The Rosseland mean opacity to absorption, in cm^2/g

The Planck, reciprocal Planck, and Rosseland mean opacities to extinction are
defined as

.. math:: \chi_{\rm P} = \frac{\int \chi_\nu B_\nu\,d\nu}{\int B_\nu\,d\nu}
          \quad\quad
          \chi_{\rm P}^{\rm inv} = \frac{\int B_\nu\,d\nu}{\int \chi_\nu^{-1} B_\nu\,d\nu}
          \quad\quad
          \chi_{\rm R} = \frac{\int \partial B_\nu / \partial T\,d\nu}{\int \chi_\nu^{-1}\,\partial B_\nu / \partial T\,d\nu}

where :math:`B_\nu` is the Planck function at the temperature :math:`T`, and
similarly for the mean opacities to absorption, with :math:`\kappa_\nu`
instead of :math:`\chi_\nu`. Note that the reciprocal Planck mean opacity is
therefore not the reciprocal of the Planck mean opacity, but the reciprocal
of the Planck mean of :math:`1/\chi_\nu`. The reciprocal Planck mean
opacity is the mean opacity that appears in the diffusion coefficient for the
modified random walk: it is defined in Equation 2 of the `corrigendum`_ to
`Robitaille (2010)`_, where it is called the Planck reciprocal mean opacity
(in Equation 18 of the original paper, it was incorrectly identified as the
Rosseland mean opacity). The Rosseland mean opacities are used for example by
the partial diffusion approximation.

.. _Robitaille (2010): https://doi.org/10.1051/0004-6361/201015025
.. _corrigendum: https://doi.org/10.1051/0004-6361/201015025e

The temperatures specified should range from 0.1K (or less) to a
temperature safely above the maximum temperature expected for the dust in
the system.

``emissivities``
----------------

This dataset should consist of a table specifying the emissivities. The
columns should be:

* ``nu``: The frequencies at which the emissivity is specified.

* ``jnu``: The emissivity for the specified frequency, as a function of
  specific energy. This should be a vector column, where the width of the
  column is the number of values specified in the ``emissivity_variable``
  dataset.

The emissivities do not need to be normalized, since Hyperion normalizes the
emissivities for each value of the specific energy when reading in the dust
file, so only the spectral shape of the emissivities matters.

``emissivity_variable``
-----------------------

This dataset should consist of a single-column table. The column should be
``specific_energy`` and should give the specific energies for which the
emissivities are tabulated in the ``emissivities`` dataset.
