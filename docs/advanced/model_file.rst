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

Each group should have the layout described in :doc:`dust_file` or should be an external HDF5 link to an actual dust file.

Grid
====

The **Grid** group should contain two sub-groups, **Geometry**, and **Physics**, which are described below.

Geometry
--------

This group describes the geometry of the model (what type of grid is used, and
the position of the cell walls). The group should have an attribute, ``grid_type``, giving the type of grid as a string which can be:

* ``car``: cartesian grid
* ``sph_pol``: spherical polar grid
* ``cyl_pol``: cylindrical polar grid
* ``amr``: adaptive mesh refinement grid (AMR)
* ``oct``: octree grid

The content of the group then depends on the type of grid:

Cartesian (``car``)
^^^^^^^^^^^^^^^^^^^

The **Geometry** group should contain three tabular datasets named ``walls_1``,
``walls_2``, and ``walls_3``, which should each contain one column. The
``walls_1`` table should contain a column ``x`` giving the x position of the
cell walls as floating point values. Similarly, ``walls_2`` and ``walls_3``
should contain one column each, named ``y`` and ``z`` respectively, giving the
y and z position of the grid cell walls.

Spherical Polar (``sph_pol``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **Geometry** group should contain three tabular datasets named ``walls_1``,
``walls_2``, and ``walls_3``, which should each contain one column. The
``walls_1`` table should contain a column ``r`` giving the radial position of
the cell walls as floating point values. Similarly, ``walls_2`` and ``walls_3``
should contain one column each, named ``t`` and ``p`` respectively, giving the
theta and phi position of the grid cell walls.

Cylindrical Polar (``cyl_pol``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **Geometry** group should contain three tabular datasets named ``walls_1``,
``walls_2``, and ``walls_3``, which should each contain one column. The
``walls_1`` table should contain a column ``w`` giving the radial position of
the cell walls as floating point values. Similarly, ``walls_2`` and ``walls_3``
should contain one column each, named ``z`` and ``p`` respectively, giving the
z and phi position of the grid cell walls.

AMR (``amr``)
^^^^^^^^^^^^^

The **Geometry** group should contain an attribute ``nlevels`` giving the
number of levels in the grid, as an integer, as well as one sub-group per
level. These groups should be formatted as **level_%05d** (i.e.
**level_00001**, **level_00002**, etc.) starting at **level_00001**.

Each **level_*** group should then contain an attribute ``ngrids`` giving the number of grids in the level, as an integer, as well as one sub-group per grid in the level. The sub-groups should be formatted as ``grid_%05d`` (e.g. **grid_00001**, **grid_00002**) starting at **grid_00001**.

Each **grid_*** group should contain the following attributes:

* ``xmin``, ``xmax``, ``ymin``, ``ymax``, ``zmin``, ``zmax``: the boundaries of
  the grid, as floating point values.

* ``n1``, ``n2``, ``n3``: the number of cells (not walls) in each direction, as
  integers.

Octree (``oct``)
^^^^^^^^^^^^^^^^

The **Geometry** group should contain the following attributes:

* ``x``, ``y``, ``z``: the coordinates of the center of the parent cell, as floating point values, in cm

* ``dx``, ``dy``, ``dz``: the size of the parent cell, as floating point values, in cm

In addition, the group should contain a 1-d array representing the ``refined`` array described in :doc:`indepth_oct`. The array should be given as an integer array instead of a boolean array.

Physics
-------

This group describes the input quantities such as the density and optionally
the specific energy of the dust. In all cases where a ``density`` array should be specified, you may also give a ``specific_energy`` array with the same dimensions - this can be used as the initial temperature, or can be used as the temperature to use for the images/SED if the number of temperature iterations is set to zero.

Cartesian
^^^^^^^^^

The **Physics** group should contain a 4-d dataset array named ``density`` giving the density in c.g.s in each cell. The dimensions of the array should be ``(n_dust, n_z, n_y, n_x)``.

Spherical Polar
^^^^^^^^^^^^^^^

The **Physics** group should contain a 4-d dataset array named ``density`` giving the density in c.g.s in each cell. The dimensions of the array should be ``(n_dust, n_p, n_t, n_r)``.

Cartesian
^^^^^^^^^

The **Physics** group should contain a 4-d dataset array named ``density`` giving the density in c.g.s in each cell. The dimensions of the array should be ``(n_dust, n_p, n_z, n_w)``.

AMR
^^^

The **Physics** group should contain a structure similar to that used to
represent the geometry. The ``nlevels`` and ``ngrids`` attributes are not
needed, only the nested **level_*** and **grid_*** groups. Each **grid_***
group should then contain a 4-d dataset array named ``density`` giving the
density in c.g.s in each cell. The dimensions of the array should be ``(n_dust,
n_z, n_y, n_x)``.

Octree
^^^^^^

The **Physics** group should contain a 1-d dataset array named ``density`` giving the density in c.g.s in each cell. Each cell in this array should match a cell in the ``refined`` array discussed in `Geometry`_.

Sources
=======

This should contain one group per source. The name of the
groups is not important, and the Python code uses names such as
**source_00001**. Each sub-group will contain certain attributes and datasets depending on the source type, as described below.

Common attributes
-----------------

All sources should have the following attributes:

* ``type``: the type of the source, given as a string. This can be ``point``
  (for point sources), ``sphere`` (for spherical sources), ``map`` (for diffuse
  sources), ``extern_sph`` (for external spherical illumination),
  ``extern_box`` (for external illumination from a box), or ``plane_parallel``
  (for a plane-parallel beam).

* ``luminosity``: the luminosity of the source, as a floating point value, in
  c.g.s

* ``peeloff``: whether to include the source when computing images with
  peeling-off, as a string that should be either ``yes`` or ``no``.

* ``spectrum``: the type of spectrum to use, as a string. This can be either:

  * ``spectrum``, to indicate that a spectrum has been numerically specified.
    In this case, the group representing the source should also contain a
    tabular dataset with two columns: ``nu``, the frequency in Hz, and ``fnu``,
    the monochromatic flux per frequency (the exact units are not important,
    because the spectrum is renormalized).

  * ``temperature``, to specify that a temperature has been specified. In this
    case, the temperature should be given as an attribute ``temperature``, as a
    floating-point value.

  * ``lte``, to indicate that the source should emit from the dust emissivities
    in the cell (this is used mainly for diffuse sources). In this case, no
    addition attributes need to be specified.

Point sources (``point``)
-------------------------

A group representing a point source should contain the following attributes in
addition to the `Common attributes`_ discussed above:

* ``x``, ``y``, and ``z``: the cartesian position of the source, as floating
  point values, in cm

Spherical sources (``sphere``)
------------------------------

A group representing a spherical source should contain the following attributes
in addition to the `Common attributes`_ discussed above:

* ``x``, ``y``, and ``z``: the cartesian position of the center of the source,
  as floating-point values, in cm

* ``r``: the radius of the sphere, as a floating point value, in cm

* ``limb``: whether to include limb darkening, as a string that can be either
  ``yes`` or ``no``.

.. TODO: mention spots

Diffuse sources (``map``)
-------------------------

In addition to the `Common attributes`_ discussed above, a group representing a diffuse source should contain a dataset called ``Luminosity map`` containing the relative luminosity of each cell as a 3-d array. The dimensions of the grid should be identical to the density grid (see `Grid`_).

External spherical sources (``extern_sph``)
-------------------------------------------

A group representing external illumination from a spherical source should
contain the following attributes in addition to the `Common attributes`_
discussed above:

* ``x``, ``y``, and ``z``: the cartesian position of the center of the source,
  as floating-point values, in cm

* ``r``: the radius of the sphere, as a floating point value, in cm

External box sources (``extern_box``)
-------------------------------------

A group representing external illumination from a box source should contain the
following attributes in addition to the `Common attributes`_ discussed above:

* ``xmin``, ``xmax``, ``ymin``, ``ymax``, ``zmin``, ``zmax``: the lower and upper bounds definining the box, as floating-point values, in cm.

Plane parallel sources (``plane_parallel``)
-------------------------------------------

A group representing a plane-parallel beam source should contain the following
attributes in addition to the `Common attributes`_ discussed above:

* ``x``, ``y``, and ``z``: the cartesian position of the center of the source,
  as floating-point values, in cm

* ``r``: the radius of the sphere, as a floating point value, in cm

* ``theta``, ``phi``: the 3-d angle giving the direction of the beam in
  spherical polar coordinates, as floating point values, in degrees.

.. TODO: mention point source collection

Output
======

The ``Output`` group should have four attributes `output_density`,
`output_density_diff`, `output_n_photons`, and `output_specific_energy`, which
should be set to a string to indicate whether to output the quantity after the
last iteration (``last``), after every iteration (``all``), or never
(``none``). The ``density_diff`` quantity is the density difference compared to
the input density (which will be non-zero in cases for example where one uses
dust sublimation).

In addition, the ``Output`` group should contain two sub-groups, ``Binned`` and
``Peeled``, that can be used to specify the parameters of the output
images/SEDs. Both groups should always be present, even if they are empty. The content of these groups is described in the following two sections:

Binned
------

If you want to compute images using binning of escaping photons (not
recommended in most cases as it is inefficient and causes angle-averaging of
outputs), then set the ``n_theta`` and ``n_phi`` parameters, which should be
used to indicate, as integers, the number of bins in the theta and phi
directions respectively.

Peeled
------

This group should contain as many sub-groups as image/SED sets you want to
compute (the name of the sub-groups is unimportant). Each sub-group should then
contain the following attributes:

* ``n_view``: the number of viewing angles for the image/SED, given as an
  integer.

* ``compute_image``: whether to compute images, given as a string that can be
  ``yes`` or ``no``. If this is ``yes``, then the following attributes should
  also be specified:

  * ``n_x`` and ``n_y``: the dimensions of the image, as integers

  * ``x_min``, ``x_max``, ``y_min``, and ``y_max``: the lower and upper bounds
    of the image as floating point values, in cm

* ``compute_sed``: whether to compute SEDs, given as a string that can be
  ``yes`` or ``no``. If this is ``yes``, then the following attributes should
  also be specified:

  * ``n_ap``: the number of apertures to compute the SED for

  * ``ap_min``, ``ap_max``: the smallest and largest apertures to use. If
    ``n_ap`` is 1, then ``ap_min`` should be the same as ``ap_max``.

* ``track_origin``: indicates whether the origin of the photon (e.g. emission
  vs scattering, or which source it originated from) should be retained in the
  output image. This can be:

  * ``no``: no photon tracking is done

  * ``basic``: photons are split into ones emitted or scattered, and whether
    they were last emitted from a source or from the dust.

  * ``detailed``: as for ``basic``, but also keeping the ID of the source or
    dust population last emitted from.

  * ``scatterings``: photons are split into ones emmitted by sources or dust,
    and split by the number of times they scatter.

* ``track_n_scat``: an integer giving the maximum number of scatterings to
  record if ``track_origin`` is set to ``scatterings``.

* ``uncertainties``: whether to compute and output uncertainties on the images
  and/or SEDs. This should be specified as a string that can be ``yes`` or
  ``no``.

* ``n_wav``: the number of wavelengths/frequencies to compute the images and/or
  SEDs for.

  * If using monochromatic radiative transfer, then the minimum and maximum
    frequency of the image and/or SED should be specified with two attributes
    ``inu_min`` and ``inu_max``, which should be given as integers giving the
    indices to the ``frequencies`` array (the first array element should be 1).

  * If not using monochromatic radiative transfer, then the minimum and maximum
    wavelengths of the image and/or SED should be specified with two attributes
    ``wav_min`` and ``wav_max``, which should be given as floating point
    values, in microns.

* ``d_min`` and ``d_max``: these give the minimum and maxium depth within which
  to use photons for the image/SED. Unless you need this and understand the
  implications, you should set this to ``-inf`` and ``+inf`` respectively if
  ``inside_observer`` is ``no``, and ``0`` and ``+inf`` respectively if
  ``inside_observer`` is ``yes``.

* ``inside_observer``: whether to compute the image from an observer located
  inside the grid, or from an observer at infinity. This should be given as a
  string that can be either ``yes`` or ``no``. In most cases you will likely
  want to use ``no``.

* ``ignore_optical_depth``: whether to ignore optical depth effects when
  computing the image/SED. This should be given as a string that can be either
  ``yes`` or ``no``, and should be set to ``no`` in most cases. This can be
  useful for debugging and for understanding how much optical depth effects are
  affecting an image or SED.

* If ``inside_observer`` is ``yes``, then the position of the observer should
  be given with the ``observer_x``, ``observer_y``, and ``observer_z``
  attributes, as floating point values, in cm.

* If ``inside_observer`` is ``no``, then the origin for the peeling-off should
  be given with the ``peeloff_x``, ``peeloff_y``, and ``peeloff_z`` attributes,
  as floating point values, in cm. In most cases, these should be set to zero.

In addition, the group should contain a table dataset with two columns,
``theta`` and ``phi``, giving the viewing angles as floating point values, in
degrees.

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


* ``forced_first_interaction``: whether to use the forced first interaction
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

