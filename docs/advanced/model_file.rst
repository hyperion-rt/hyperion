=======================
Model Input HDF5 Format
=======================

.. _input_intro:

Introduction
============

The radiation transfer code reads in an HDF5 file as input. The contents of the file have to follow a specific layout. To differentiate this format from other HDF5 files, we use the ``.rtin`` extension for input files to the radiation transfer code.

An ``.rtin`` file should contain the following four groups::

    Dust/
    Grid/
    Sources/
    Output/

These are described in `Dust`_, `Grid`_, `Sources`_, and `Output`_
respectively. In addition, a number of attributes should be set at the root level, and these are described in `Root-level attributes`_.

Dust
====

The `Dust`_ group should contain as many groups (or external/internal links
to groups) as dust types. The groups should be named::

    dust_001/
    dust_002/
    ...

Each group should have the layout described in :doc:`dust_file`.

Grid
====

This section is incomplete.

Sources
=======

This section is incomplete.

Output
======

This section is incomplete.

Root-level attributes
=====================

The overall configuration for the model should be specified as attributes for the root group in the HDF5 file. The parameters needed are described in the following sub-sections.

General
-------

* ``python_version``: the version of the Hyperion Python library used to
  generate the file. If you are not generating the file with the Hyperion
  Python library (which is probably the case if you are reading this page)
  then set it to '0.8.7' since that is the version for which the format in
  this page is described.

* ``physics_io_bytes``: whether to write out the physical quantities using
  4- or 8-byte floating point values. Should be either ``4`` or ``8``
  (integer).

Iterations
----------

* ``n_lucy_iter``: Number of temperature-calculating Lucy iterations
  (integer)

* ``check_convergence``: Whether to check for convergence in the specific
  energy calculation. Should be ``yes`` or ``no`` (string).

* ``convergence_absolute``: the threshold for absolute changes in the
  specific energy (float).

* ``convergence_relative``: the threshold for relative changes in the
  specific energy (float).

* ``convergence_percentile``: the percentile at which to check the absolute
  and relative changes in the specific energy (float).

See :ref:`convergence` for the latter three.

Diffusion approximation
-----------------------

* ``mrw``: Whether or not to use the modified random walk (MRW) algorithm.
  Should be ``yes`` or ``no`` (string).

* ``pda``: Whether or not to use the partial diffusion approximation (PDA)
  algorithm. Should be ``yes`` or ``no`` (string).

If ``mrw`` is ``yes``, the following two attributes should be set:

* ``mrw_gamma``: The gamma parameter for the modified random walk as
  described in :ref:`diffusion` (float).

* ``n_inter_mrw_max``: The maximum number of MRW interactions before a
  photon is killed (integer).

Images/SEDs
-----------

* ``raytracing``: Whether to do a raytracing iteration at the end of
  the calculation. Should be ``yes`` or ``no`` (string).

* ``monochromatic``: Whether to calculate the images/SEDs in monochromatic
  mode. Should be ``yes`` or ``no`` (string).

Number of photons
-----------------

The following attributes are required:

* ``n_stats``: how often to display performance stats. For the MPI-enabled
  code, this also determines the chunk of photons to dispatch to each thread
  at a time (integer).

If ``n_initial_iter`` is non-zero, then the following photon number should be specified:

* ``n_initial_photons``: number of photons to emit per iteration in the
  initial iterations (integer).

If ``raytracing`` is ``yes``, then the following photon numbers should be specified:

* ``n_ray_photons_sources``: number of raytracing photons from sources
  (integer). Does not need to be specified if there are no sources.

* ``n_ray_photons_dust``: number of raytracing photons from dust (integer).
  Does not need to be specified if there are no dust density grids.

If ``monochromatic`` is ``yes``, then the following photon numbers should be specified:

* ``n_last_photons_sources``: the number of photons (per frequency) to emit
  from sources in the imaging iteration (integer). Does not need to be
  specified if there are no sources.

* ``n_last_photons_dust``: the number of photons (per frequency) to emit
  from dust in the imaging iteration (integer). Does not need to be
  specified if there are no dust density grids.

Miscellaneous
-------------


* ``forced_first_scattering``: whether to use the forced first scattering
  algorithm. Should be one of ``yes`` or ``no`` (string).

* ``kill_on_absorb``: whether to kill photons when they are absorbed rather
  than re-emitting them (useful for scattering-only calculations). Should be
  one of ``yes`` or ``no`` (string).

* ``n_inter_max``: the maximum number of interactions a photon can have
  before being it is killed (integer).

* ``n_reabs_max``: the maximum number of times a photon can be re-absorbed
  before it is killed (integer).

Optional
--------

The following attributes are optional:

* ``sample_sources_evenly``: whether to emit the same number of photons from
  each source (as opposed to emitting a number of photons proportional to
  the luminosity). Should be ``yes`` or ``no`` (string). Defaults to ``no``.

* ``enforce_energy_range``: whether to always reset values below the minimum
  and above the maximum specific energy to the bounds of the range. Should
  be ``yes`` or ``no`` (string). Defaults to ``yes``.

