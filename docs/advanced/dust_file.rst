=================
Dust HDF5 Format
=================

Overview
========

Hyperion requires information on the dust properties, such as the albedo,
opacity, mean opacities, and emissivities. These need to be packaged in an
HDF5 format that is described in the :ref:`specification`. However, in most
cases you do not need to create these files from scratch.

Creating dust files
===================

Converting dust files from the ``ttsre`` code
----------------------------------------------

Converting a dust file from the format used by Barb Whitney's ``ttsre`` code
is straightforward. Given for example the ``kmh.par`` file, you can create an
HDF5 file in the format required for Hyperion using::

    from hyperion.dust import SimpleSphericalDust
    d = SimpleSphericalDust('kmh.par')
    d.write('kmh.hdf5')

This is a one-time operation for each dust type - once the HDF5 file has been
created, you do not need to recreate it every time you want to set up a model.
You will only need to run this again if explicitly asked to after a Hyperion
update.

.. note:: If you get an error that looks like::

              Exception: x values are out of interpolation bounds

          this means that the optical properties are not defined for a wide
          enough range of wavelengths. You can extrapolate the existing values
          using::

              d.optical_properties._extrapolate(1.e-3, 1.e7)

          before the ``d.write`` command.

In addition, you can plot an overview of the optical properties using::

    d.plot('kmh.png')

.. _specification:

Dust file HDF5 format specification
===================================

An HDF5 dust file should contain 5 datasets. The root of the file should contain the following attributes:

* ``emissvar``: whether the emissivity is specified as a function of specific
  energy absorbed in each cell (``E``) or another quantity (but this is not
  supported at this time).

* ``version``: this should be set to ``1`` - the version described in this
  section.

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

* ``specific energy``: The specific energy for which the mean opacities are
  given

* ``chi_planck``: The Plank mean opacity to extinction, in cm^2/g

* ``chi_rosseland``: The Rosseland mean opacity to extinction, in cm^2/g

* ``kappa_planck``: The Plank mean opacity to absorption, in cm^2/g

* ``kappa_rosseland``: The Rosseland mean opacity to absorption, in cm^2/g

The temperatures specified should range from 0.1K (or less) to a
temperature safely above the maximum temperature expected for the dust in
the system.

``emissivities``
----------------

This dataset should consist of a table specifying the emissivities. The
columns should be:

* ``nu``: The frequencies at which the emissivity is specified.

* ``j_nu``: The emissivity for the specified frequency, as a function of
  specific energy. This should be a vector column, where the width of the
  column is the number of values specified in the ``emissivity_variable``
  dataset.
  
``emissivity_variable``
-----------------------

This dataset should consist of a single-column table. The column should be
``specific_energy`` and should give the specific energies for which the
emissivities are tabulated in the ``emissivities`` dataset.
