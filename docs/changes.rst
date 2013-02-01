Version History
===============

0.9.2 (unreleased)
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

Improvements
^^^^^^^^^^^^

- PyFITS, PyWCS, and ATpy are no longer required for Hyperion. Instead, the
  `Astropy <http://www.astropy.org>`_ package is now required as a dependency.

- Updated download link for MPICH2

Bug fixes
^^^^^^^^^

- Fix compatibility with Numpy 1.8.0.dev

- Fix coverage testing for Python 3

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
