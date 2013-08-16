Version History
===============

0.9.2 (2013-08-16)
------------------

New Features
^^^^^^^^^^^^

- :meth:`~hyperion.model.ModelOutput.get_sed` and
  :meth:`~hyperion.model.ModelOutput.get_image` now return SED and
  Image objects that contain meta-data in addition to the data itself. For
  example, images contain information about the field of view (in
  physical/angular units, where appropriate), and information about the units
  is also included. The old syntax of ``wav, nufnu = m.get_sed(...)`` will
  still work, but the meta-data will not be accessible in those cases.

- New library of dust models, accessible in :doc:`dust/dust`

- It is now possible to read in previous models completely, including the
  density structure, geometry, sources, dust, and configuration, using the
  :meth:`~hyperion.model.Model.read` method. In addition, new methods
  individual methods :meth:`~hyperion.model.Model.use_sources`,
  :meth:`~hyperion.model.Model.use_image_config`,
  :meth:`~hyperion.model.Model.use_run_config`, and
  :meth:`~hyperion.model.Model.use_output_config` allow more detailed control
  over reading in parameters from previous models.

- It is now possible to force overwrite Hyperion output from the command-line
  using the ``-f`` option::

    hyperion -f input output

  or when using the individual fortran binaries::

    mpirun -n 8 hyperion_car_mpi -f input output

  This will likely be useful for users of computer clusters who don't want a
  job to fail just because the output file already exists.

- Regular Cartesian grids can now also be exported for viewing in `yt
  <http://yt-project.org/>`_ (as was previously possible for AMR and Octree
  grids).

- A new function, :func:`~hyperion.model.helpers.run_with_vertical_hseq`,
  is available to help with the calculation of vertical Hydrostatic
  equilibrium in disks. Note that this feature is still experimental and
  should be used with care.

- A new function, :func:`~hyperion.model.helpers.tau_to_radius`, is available
  to compute, for spherical polar grids, the optical depth from infinity to a
  given radius.

Improvements
^^^^^^^^^^^^

- PyFITS, PyWCS, and ATpy are no longer required for Hyperion. Instead, the
  `Astropy <http://www.astropy.org>`_ package is now required as a dependency.

- Updated download link for MPICH2

- The ``rho_0`` attribute for disks is now a property, not a method, and can
  be set by the user instead of the disk mass.

- The documentation has been improved and fixed in places thanks to user
  feedback.

- AnalyticalYSOModel instances are no longer 'static' once they have been
  written out (this means one can write out a model, change a parameter, and
  write out a new different model, which was not possible previously).

- The Fortran code now reads in dust models faster because it computes all
  cumulative distribution functions more efficiently.

- Statistics for killed photons are now kept for each iteration rather than
  just summing all of them.

Bug fixes
^^^^^^^^^

- Fix compatibility with Numpy 1.8.0.dev

- Fix coverage testing for Python 3

- Fixed an issue which caused temporary files to not be deleted after running
  tests.

API changes
^^^^^^^^^^^

- The ``AnalyticalYSOModel.evaluate_optically_thin_radii()`` method has been
  removed.

0.9.1 (2012-10-26)
------------------

New Features
^^^^^^^^^^^^

- Updated hyperion2fits to extract binned images

- Added wmax= option for AnalyticalYSOModel.set_cylindrical_grid_auto

Improvements
^^^^^^^^^^^^

- Made deps/fortran/install.py script more robust to architecture, and to lack
  of zlib library.

- Ensure that spectrum always gets converted to floating-point values

- Give a more explicit error message if optical properties for dust are not
  set.

Bug fixes
^^^^^^^^^

- Fixed bug that prevented BipolarCavity from being used

- Ensure that get_quantities works even if no initial iterations were computed

- Fix scattering for cases where P2=0. The code could sometimes crash if a mix
  of isotropic and non-isotropic dust was used (reported by M. Wolff).

- Fix a bug that occurred when outputting multiple images with the depth
  option (reported and fixed by T. Bowers) [#21, #22]

0.9.0 (2012-07-27)
------------------

- Initial public release.
