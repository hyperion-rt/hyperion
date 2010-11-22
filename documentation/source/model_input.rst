=======================
Model Input HDF5 Format
=======================

.. _input_intro:

Introduction
============

The radiation transfer code reads in an HDF5 file as input. The contents of the file have to follow a specific layout. To differentiate this format from other HDF5 files, we use the ``.rtin`` extension for input files to the radiation transfer code.

An ``.rtin`` file should follow this layout::

    model.rtin/
        Dust/
            Dust types
            [plus one group per dust type]
            ...
        Grid/
            Geometry
                Walls 1
                Walls 2
                Walls 3
            Physics
                Density
                Temperature [optional]
        Output/
            Binned
                Group 00001 [optional]
            Peeled
                Group 00001 [optional]
                Group 00002 [optional]
                ...
        Sources/
            Source 1
            Source 2
            ...

The following sections describe in detail each of these groups.

Overall configuration
======================

The overall configuration for the model should be specified as attributes for the root group in the HDF5 file. The parameters needed are described in the following sub-sections.

Iterations
----------

* ``n_lucy_iter``: Number of temperature-calculating Lucy iterations.

* ``diffusion``: Whether or not to use the diffusion approximation in
  optically thick regions. Should be ``yes`` or ``no``.

* ``raytracing``: Whether or not to do a raytracing iteration at the end of
  the calculation. Should be ``yes`` or ``no``.

Number of photons
-----------------

The following parameters control the number of photons for various purposes

* ``n_stats``: how often to display performance stats. For the MPI-enabled
  code, this also determines the chunk of photons to dispatch to each thread
  at a time.

* ``n_lucy_photons``: The number of photons to use during each iteration in
  the iterations used to calculate the temperature using the Lucy algorithm.

* ``n_last_photons``: The number of photons to use during the iteration
  following the Lucy iterations, when images/SEDs are calculated.

* ``n_ray_photons``: The number of photons to use during the extra raytracing
  iteration.

Output
------

* ``make_binned_images``: Whether to compute binned images/SEDs. Should be
  ``yes`` or ``no``.

* ``binned_config``: The name of the file inside ``output/`` that contains the
  parameters for the binned images/SEDs.

* ``make_peeled_images``: Whether to compute peeloff images/SEDs. Should be
  ``yes`` or ``no``.

* ``n_peeled``: Number of groups of peeloff images (one group can contain
  multiple viewing angles).

* ``peeled_config(?)``: The name of the file inside ``output/`` that contains
  the parameters for the peeloff images/SEDs. The question mark should be
  replaced by the group number. There should be ``n_peeled`` times
  ``peeled_config`` entries.

Dust parameters
===============

Dust list
---------

The ``dust/`` directory should contain a file named ``dust.txt`` that lists the dust to use for each dust index, and the dust files themselves. The ``dust.txt`` file should contain two columns::

    001 www003.fits          
    002 kmh.fits             
    003 r400_ice095.fits     
    004 kmh.fits             
    ...
    
The first column is the dust index, and should be in the correct order. This is redundant because it is equal to the line number, but is included for safety. The second column is the filename of the dust file. All files listed here should be included in the ``dust/`` directory, and should be in FITS format with the layout specified below.

Grid parameters
===============

The ``grid/`` directory should contain two files, one describing the geometry of the grid, and one specifying the densities in the grid. Both files should be in the FITS format.

Geometry file (3D cartesian grid)
---------------------------------

The FITS geometry file should contain 6 HDUs (including the primary HDU). If possible, all values should be specified in 8-byte precision. In the following description, we use ``nx``, ``ny`` and ``nz`` to refer to the number of **grid cells** in the x, y, and z directions.

HDU 0 (primary)
^^^^^^^^^^^^^^^

This HDU should contain a 3D array with dimensions ``nx x ny x nz``, containing the volumes of each cell in cm^3.

HDU 1
^^^^^

This HDU should contain a 4D array with dimensions ``nx x ny x nz x 6``, containing the areas of the 6 cell walls in cm^2. The order of the walls is xmin, xmax, ymin, ymax, zmin, and zmax.

HDU 2
^^^^^

This HDU should contain a 4D array with dimensions ``nx x ny x nz x 3``, containing the width of the cell in the 3 directions in cm. The order of the widths is dx, dy, and dz.

HDU 3
^^^^^

This HDU should contain a single-column binary table. The column should be ``x`` and give the ``nx + 1`` wall positions in the x-direction.

HDU 4
^^^^^

This HDU should contain a single-column binary table. The column should be ``y`` and give the ``ny + 1`` wall positions in the y-direction.

HDU 5
^^^^^

This HDU should contain a single-column binary table. The column should be ``z`` and give the ``nz + 1`` wall positions in the z-direction.

Geometry file (3D spherical polar grid)
---------------------------------------

The 3D spherical polar grid geometry file is identical to the cartesian grid, where ``x`` should be replaced by ``r``, ``y`` by ``t`` (theta), and ``z`` by ``p`` (phi).

Geometry file (adaptive grid)
-----------------------------

Not implemented

Source parameters
=================

The ``sources/`` directory should contain a file named ``sources.txt`` that lists all the sources of photons, and any ancillary files. The ``sources.txt`` file should contain one line per source. All lines should start with a string giving the type of the source, followed by the luminosity of the source, in ergs/cm^2/s. The format of the rest of the line then depends on the source type, as described below.

Source list
-----------

Point Source (``point``)
^^^^^^^^^^^^^^^^^^^^^^^^

The columns should be:

* Column ``3`` to ``5``: The x, y, and z position of the source (in cm)

* Column ``6``: The filename of the spectrum to use, prefixed by ``spec:``
  (e.g. ``spec:xray.fits``), or the temperature of a blackbody, prefixed by
  ``temp:`` (e.g. ``temp:10000``).

Spherical Source (``sphere``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The columns should be:

* Columns ``3`` to ``5``: The x, y, and z position of the source (in cm)

* Column ``6``: The radius of the source (in cm)

* Column ``7``: A single character indicating whether to include (``T``) or
  not (``F``) limb darkening.

* Column ``8``: The filename of the spectrum to use, prefixed by ``spec:``
  (e.g. ``spec:xray.fits``), or the temperature of a blackbody, prefixed by
  ``temp:`` (e.g. ``temp:10000``).
  
.. note::
    It is possible to add spots to a spherical source. These should be listed
    on the lines immediately following their parent spherical source, and
    should start with ``+spot``. As for normal sources, the second column
    should give the luminosity. The following columns should list the
    longitude, latitude, and size of the spot (in degrees), and the final
    column should list the spectrum or temperature using the same ``spec:`` or
    ``temp:`` syntax as for point and spherical sources.

Spectrum file format
--------------------

A FITS atmosphere file should contain 2 HDUs (including the primary HDU). If possible, all values should be specified in 8-byte precision. The HDUs should be:

HDU 0 (primary)
^^^^^^^^^^^^^^^

The primary HDU should be left empty.

HDU 1
^^^^^

This HDU should contain a binary table specifying the spectrum. The columns should be:

* ``nu``: The frequency (in Hz)

* ``wav``: The wavelength (in microns)

* ``fnu``: The flux, given as :math:`F_\nu`. The units are not important since
  the spectrum is used only as a probability distribution function. The
  luminosity is specified in the ``sources.txt`` file described above.

Output parameters
=================

The ``output/`` directory should contain configuration files (see
:ref:`input_intro` for information of the ``.conf`` format), each one
describing an image/SED group listed in the ``run.conf`` file. The required
parameters are:

* ``n_x`` and ``n_y``: Number of pixels in the binned images.

* ``x_min``, ``x_max``, ``y_min``, and ``y_max``: The range of x and y values
  for the binned images (in cm).

* ``n_wav``: Number of wavelengths in the binned images

* ``wav_min``, and ``wav_max``: The range of wavelengths values for the binned
  images (in microns).

* ``n_ap``: Number of apertures for the binned SEDs.

* ``ap_min``, and ``ap_max``: The range of aperture values for the binned SEDs
  (in cm).

Binned images/SEDs
------------------

For a configuration file used for binned images/SEDs, the ``n_theta`` and
``n_phi`` parameters is also required. These give the number of
theta/phi angle bins for the binned images/SEDs
  
Peeling-off images/SEDs
-----------------------

For a configuration file used for peeling-off images/SEDs, the ``n_view`` parameter is also required. This gives the number of viewing angles to compute images/SEDs for. Finally, the angles should be specified using the ``theta`` and ``phi`` angles, e.g::

    theta = 10. 20. 30. 40.
    phi   =  0.  0.  0.  0.
    
The number of viewing angles specified should match ``n_view``


