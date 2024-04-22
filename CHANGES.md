## v0.9.10 - 2019-02-08

- Fixed compatibility with HDF5 1.10.x. [#188]
- Simplify installation instructions, remove custom installer for Fortran dependencies, and remove CMake build system. [#206]
- Fixed compatibility with latest version of h5py. [#201]
- Fixed bug in method to truncate scattering matrix in dust properties. [#199]

## v0.9.9 - 2017-08-28

- Added a workaround for a bug in ifort v16 that prevented compilation. [#190]
- Improved the accuracy of the Wood and Reynolds (1999) forced first scattering algorithm for very small optical depths, and implemented the Baes et al - 2016) composite biasing algorithm for forced first interaction. In addition, the `Model.set_forced_first_scattering` method has been renamed to `Model.set_forced_first_interaction`, and this method now optionally takes the algorithm to use (`wr99` or `baes16`) and in the case of the `baes16` algorithm, `baes16_xi` can be passed. [#189]
- Improved the efficiency of monochromatic radiative transfer: instead of propagating photons until they are absorbed, photon packets are forced to scatter and their energy is lowered at every interaction by the albedo until the energy goes below a threshold. This threshold defaults to 1.e-10, but this can be changed using the `energy_threshold` argument to `Model.set_monochromatic`. [#158]
- Remove bundled numpydoc and fix Python 3.6 compatibility. [#192]

## v0.9.8 - 2016-09-27

- Refactored installation instructions to provide a more streamlined experience. [#161]
- Fixed a severe bug in the set-up of AMR grids that caused photons to be terminated before escaping from grids in specific cases. [#167]
- Fixed bug in ordering of axes for Cartesian and AMR grids for yt 3.x. [#168]
- Add installation instructions for FreeBSD (thanks to Aaron Lauterer). [#173]
- Fixed a bug that caused a VoronoiGrid initialized from another VoronoiGrid to not include the sparse neighbor information. [#177]
- Check for NaNs when writing out dust or model HDF5 file. [#179]
- Fix the installation script for dependencies to work with GCC 5. [#182]
- Remove deprecated frequencies= argument from set_monochromatic (instead, the wavelengths= argument should be used to give the wavelengths in microns). [#183]
- The set_wavelength_index_range should now take zero-based indices, consistent with the Python/Numpy convention. Note that this will break compatibility with previous scripts, so you will need to be sure to update calls to this method if you were using it. [#183]
- Update bundled decorator.py module. [#184]
- Disabled computation of SEDs for inside observers, since the apertures are not truly circular. Users should instead set up images and then do photometry on them as appropriate. [#185]
- The code base is now Python 3-compatible, and we no longer rely on 2to3. [#185]
- Update installation instructions to mention conda as the preferred installation method. [#186]
- Add workaround for a bug in HDF5 1.8 that caused issues with models that made use of filters. [#187]

## v0.9.7 - 2015-08-22

### New features

- Added support for on-the-fly filter convolution. [#119]
- Power-law and Ulrich envelopes can now be used on cylindrical polar grids. [#136]
- Provide a way to sample random positions inside Voronoi cells. [#134, #151]
- Added the ability to load AMR grids from yt datasets. [#148]

### Bug fixes

- Correctly raise an error if photons are emitted outside the Voronoi grid. [#127]
- Avoid issues with number of photons when reading in models. [#145]
- Fix a bug that caused sublimation parameters to not be correctly read in. [#133]
- Fixed Fortran dependencies install script.

### API changes

- `set_aperture_range` has been renamed to `set_aperture_radii` to emphasize that these are radii. [#132]

### Other

- Internal refactoring of how spectra are gridded for raytracing to make it easier to implement Doppler-shifted spectra in future. [#126]
- Optimizations to memory and disk usage for Voronoi grids. [#128, #154, #156]
- Improvements to warnings and documentation. [#54, #98, #125]
- Removed deprecated mctherm2hyperion script. [#144]

## v0.9.6 - 2015-02-27

### Bug fixes

- Fixed backward compatibility with files that don't include d_min/d_max

## v0.9.5 - 2015-02-17

### New features

- Added an importer function, `construct_octree`, to convert a list of SPH particles into an Octree grid. [#67]
- Addded support for Voronoi grids. Voronoi grids can be set via the `~hyperion.model.Model.set_voronoi_grid` method, passing in the `x`, `y` and `z` coordinates of the sites as arguments. [#92]
- Added the ability to make images that split up photons as a function of how many times they were scattered (using `set_track_origin('scatterings')`). [#99]
- Added the ability to kill photons when scattered using `set_kill_on_scatter`, similarly to the existing `set_kill_on_absorb`. [#101]
- Added the ability to use input files in `ModelOutput.get_quantities`. [#106]
- Added the ability to control the interpretation of the specific energy passed to `add_density_grid` using the `set_specific_energy_type` argument. [#117]
- Added the ability to use a grid from an HDF5 file on disk without reading it into the Python code, using `use_grid_from_file`. [#116]
- Added support for cmake. [#112]
- `Image` and `SED` objects now include attributes `d_min` and `d_max` that indicate the depth of the region used to contruct the image or SED. [#121]
- Fixed a bug in the computation of the Rosseland mean opacity (it was in fact the Planck reciprocal mean opacity). Dust files have now been updated to version 2 to include both the Rosseland mean opacity and the Planck reciprocal mean opacity. [#123]

### Bug fixes

- Fixed a minor bug in the logic for killing photons that have had too many interactions. [#100]
- Fixed a bug that meant that BipolarCavity instances could not be subtracted from AmbientMedium instances. [#106]

### Other improvements

- The `to_yt()` methods are now compatible with yt 3.x (3.0.1 and later recommended). [#113]
- The `uncertainties=True` mode for `get_sed` and `get_image` has now been properly vectorized and should be faster by a factor of a few when requesting polarization results. [#114]

## v0.9.4 - 2014-01-29

### New features

- Image and SED groups now have a `set_stokes` option that allows users to specify whether to save Stokes componenets other than I. Prior to this version, all Stokes components were always saved, but this resulted in an unecessarily high memory usage in many cases, so the default is now set to `False`, and users have to explicitly set `set_stokes(True)` in order to save all Stokes components. [#61]
- It is now possible to turn off the warnings that occur when photons are killed due to too many interactions, using the `warn=True/False` option for the :meth:`~hyperion.model.Model.set_max_interactions` method (and other similar methods). [#68]

### Bug fixes

- Fix Fortran dependency installer for gfortran 4.5 and earlier
- Fixed a bug that caused models using the monochromatic radiative transfer settings to not be read in correctly by :meth:`~hyperion.model.Model.read`. [#78]

### API Changes

- When using the monochromatic radiative transfer mode, users should now use the :meth:`~hyperion.conf.PeeledImageConf.set_wavelength_index_range` method instead of :meth:`~hyperion.conf.PeeledImageConf.set_wavelength_range`. [#78]

## v0.9.3 - 2013-11-14

### New features

- For models that require many point sources with a common spectrum, a new source type (point source collection) is now available. To add a point source collection, use `source = m.add_point_source_collection()`. The `source.luminosity` attribute should be set to an array with as many elements as sources, and the `source.position` attribute should be set to a 2-d array where the first dimension matches `source.luminosity`, and with 3 elements in the second dimension (x, y, and z).
- Sources can now be given names as strings, which can then be used as an argument to `source_id` in :meth:`~hyperion.model.ModelOutput.get_sed` and :meth:`~hyperion.model.ModelOutput.get_image` (when using photon tracking).
- Improved documentation to explain better in which cases dust and total densities should be used. This is summarized in :doc:`important/important`.
- Added an option to specify the minimum (relative) radial cell spacing for the :class:`~hyperion.model.AnalyticalYSOModel` class.
- Fixed bug that prevented users from setting the grid manually with the :class:`~hyperion.model.AnalyticalYSOModel` class.
- It is now possible to include multiple ambient mediums with different dust properties (this was limited to a single ambient medium property previously).
- The :meth:`~hyperion.model.Model.add_density_grid` method can now be called with a grid view for all grid types (previously this was only possible for AMR grids).
- Added dust classes to the API documentation.
- Fixed a typo in the equation for the :class:`~hyperion.densities.AlphaDisk` class, and added definitions of the scaleheight for :class:`~hyperion.densities.AlphaDisk` and :class:`~hyperion.densities.FlaredDisk`.
- Improve the reliability of the configure script.

## v0.9.2 - 2013-08-16

### New Features

- :meth:`~hyperion.model.ModelOutput.get_sed` and :meth:`~hyperion.model.ModelOutput.get_image` now return SED and Image objects that contain meta-data in addition to the data itself. For example, images contain information about the field of view (in physical/angular units, where appropriate), and information about the units is also included. The old syntax of `wav, nufnu = m.get_sed(...)` will still work, but the meta-data will not be accessible in those cases.
- New library of dust models, accessible in :doc:`dust/dust`
- It is now possible to read in previous models completely, including the density structure, geometry, sources, dust, and configuration, using the :meth:`~hyperion.model.Model.read` method. In addition, new methods individual methods :meth:`~hyperion.model.Model.use_sources`, :meth:`~hyperion.model.Model.use_image_config`, :meth:`~hyperion.model.Model.use_run_config`, and :meth:`~hyperion.model.Model.use_output_config` allow more detailed control over reading in parameters from previous models.
- It is now possible to force overwrite Hyperion output from the command-line using the `-f` option: `hyperion -f input output`
- or when using the individual fortran binaries: `mpirun -n 8 hyperion_car_mpi -f input output`. This will likely be useful for users of computer clusters who don't want a job to fail just because the output file already exists.
- Regular Cartesian grids can now also be exported for viewing in `yt <http://yt-project.org/>`_ (as was previously possible for AMR and Octree grids).
- A new function, :func:`~hyperion.model.helpers.run_with_vertical_hseq`, is available to help with the calculation of vertical Hydrostatic equilibrium in disks. Note that this feature is still experimental and should be used with care.
- A new function, :func:`~hyperion.model.helpers.tau_to_radius`, is available to compute, for spherical polar grids, the optical depth from infinity to a given radius.

### Improvements

- PyFITS, PyWCS, and ATpy are no longer required for Hyperion. Instead, the `Astropy <http://www.astropy.org>`_ package is now required as a dependency.
- Updated download link for MPICH2
- The `rho_0` attribute for disks is now a property, not a method, and can be set by the user instead of the disk mass.
- The documentation has been improved and fixed in places thanks to user feedback.
- AnalyticalYSOModel instances are no longer 'static' once they have been written out (this means one can write out a model, change a parameter, and write out a new different model, which was not possible previously).
- The Fortran code now reads in dust models faster because it computes all cumulative distribution functions more efficiently.
- Statistics for killed photons are now kept for each iteration rather than just summing all of them.

### Bug fixes

- Fix compatibility with Numpy 1.8.0.dev
- Fix coverage testing for Python 3
- Fixed an issue which caused temporary files to not be deleted after running tests.

### API changes

- The `AnalyticalYSOModel.evaluate_optically_thin_radii()` method has been removed.

## v0.9.1 - 2012-10-26

### New Features

- Updated hyperion2fits to extract binned images
- Added wmax= option for AnalyticalYSOModel.set_cylindrical_grid_auto

### Improvements

- Made deps/fortran/install.py script more robust to architecture, and to lack of zlib library.
- Ensure that spectrum always gets converted to floating-point values
- Give a more explicit error message if optical properties for dust are not set.

### Bug fixes

- Fixed bug that prevented BipolarCavity from being used
- Ensure that get_quantities works even if no initial iterations were computed
- Fix scattering for cases where P2=0. The code could sometimes crash if a mix of isotropic and non-isotropic dust was used (reported by M. Wolff).
- Fix a bug that occurred when outputting multiple images with the depth option (reported and fixed by T. Bowers) [#21, #22]

## v0.9.0 - 2012-07-27

- Initial public release.
