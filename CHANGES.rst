0.9.1 (unreleased)
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
