.. _model:

================
Arbitrary Models
================

.. note:: The current document only shows example use of some methods, and
          does not discuss all the options available. To see these, do not
          hesitate to use the ``help`` command, for example ``help
          m.write`` will return more detailed instructions on using the
          ``write`` method.

To create a general model, you will need to first import the Model class
from the Python Hyperion module::

    from hyperion.model import Model

it is then easy to set up a generic model using::

    m = Model()

The model can then be set up using methods of the ``Model`` instance. These
are described in the following sections.

Once the model is set up, the user can write it out to the disk for use
with the Fortran radiation transfer code::

    m.write('example.rtin')

.. _grid:

Grid
====

The code currently supports five types of 3-d grids:

* Cartesian grids
* Spherical polar grids
* Cylindrical polar grids
* AMR grids
* Octree grids

The following sections show how the different kinds of grids should be set up.

Regular 3-d grids
----------------

In the case of the cartesian and polar grids, the user should define the wall
position in each of the three directions, using cgs units for the spatial
coordinates, and radians for the angular coordinates. These wall positions
should be stored in one 1-d NumPy array for each dimension, with one element
more than the number of cells defined. The walls can then be used to create a
coordinate grid using methods of the form ``set_x_grid(walls_1, walls_2,
walls_3)``. The following examples demonstrate how to do this for the various
grid types

* A 10x10x10 cartesian grid from -1 pc to +1 pc in each direction::

    x = np.linspace(-pc, pc, 11)
    y = np.linspace(-pc, pc, 11)
    z = np.linspace(-pc, pc, 11)
    m.set_cartesian_grid(x, y, z)

* A 2-d 399x199 spherical polar grid::

    r = np.logspace(np.log10(rsun), np.log10(100*au), 400)
    theta = np.linspace(0., pi., 199)
    phi = np.array([0., 2*pi])
    m.set_spherical_polar_grid(r, theta, phi)

* A 3-d 100x100x10 cylindrical polar grid::

    w = np.logpsace(np.log10(rsun), np.log10(100*au), 101)
    z = np.linspace(-10*au, 10*au, 101)
    phi = np.linspace(0, 2*pi, 11)
    m.set_cylindrical_polar_grid(w, z, phi)

AMR grids
---------

AMR grids have to be constructed using the ``AMRGrid`` class::

    from hyperion.grid import AMRGrid
    amr = AMRGrid()

Levels can be added with::

    level = amr.add_level()

And grids can be added to a level with::

    grid = level.add_grid()

Grid objects have the following attributes which should be set:

* ``xmin`` - lower x position of the grid
* ``xmax`` - upper x position of the grid
* ``ymin`` - lower y position of the grid
* ``ymax`` - upper y position of the grid
* ``zmin`` - lower z position of the grid
* ``zmax`` - upper z position of the grid
* ``nx`` - number of cells in x direction
* ``ny`` - number of cells in y direction
* ``nz`` - number of cells in z direction

Once we have an AMR grid object, which we call ``amr`` here, the geometry can be set using::

    m.set_amr_grid(amr)

The quantity contained in the grid is unimportant for this step, as long as the geometry is correct.

For more details on how to create or read in an AMR object, see :ref:`amr_indepth`.

Octree grids
------------

An `Octree <http://en.wikipedia.org/wiki/Octree>`_ is a hierarchical grid format where each cell can be divided into eight children cells. At the top level is a single cell that covers the whole spatial domain being considered. To set up an Octree, the following information is needed:

* ``x``, ``y``, ``z`` - the coordinates of the center of the parent cell
* ``dx``, ``dy``, ``dz`` - the size of the parent cell
* ``refined`` a 1-d sequence of booleans giving the structure of the grid.

The ``refined`` sequence contains all the information regarding the hierarchy of the grid, and is described in :ref:`indepth_oct`. Once this sequence is set, the geometry can be set with::

    m.set_octree_grid(x, y, z, dx, dy, dz, refined)

Density and Specific Energy
===========================

Once a regular grid is set up, it is straightforward to add one or more density grids. In this step, a dust file in HDF5 format is also required. See :ref:`dustfile` for more details about creating and using dust files in HDF5
format.

Regular 3-d grids
----------------

For regular cartesian and polar grids, a 3-d NumPy array containing
the density array is required. A density grid is added with::

    m.add_density_grid(density_array, dust_file)

For example::

    m.add_density_grid(np.ones(100,100,10), 'kmh.hdf5')

This command can be called multiple times if multiple density arrays are
needed (for example if different dust sizes have different spatial
distributions).

Optionally, a specific energy distribution can also be specified using a 3-d NumPy
array using the ``specific_energy=`` argument::

    m.add_density_grid(density_array, dust_file, specific_energy=specific_energy_array)

.. note:: Specifying a specific energy distribution is only useful if the
          number of initial iterations for the RT code is set to zero (see
          `Specific Energy Calculation`_), otherwise the input specific energy
          will be overwritten with the self-consistently computed one.

AMR grids
---------

Since AMR grids have a more complex structure than regular 3-d arrays, the density should be added using an AMR object (as described in :ref:`grid`), with the quantity being used specified as an item::

    m.add_density_grid(amr_object['density'], dust_file)

for example::

    m.add_density_grid(amr['density'], 'kmh.hdf5')

For the above to work, the density should be set for each grid in each level, in a ``quantities`` dictionary that has a ``'density'`` key, which in turns contains a 3-d array with the density (see :ref:`amr_indepth` for more details).

Specific energies can be specified using the same kinds of objects and using the `specific_energy` argument::

    m.add_density_grid(amr['density], dust_file,
                       specific_energy=amr['specific_energy'])

Octree grids
------------

In the case of Octrees, densities (and optionally specific energies) should be specified in the same manner as the regular grids, but should be specified as a 1-d Numpy array with the same length as the ``refined`` list, where each density corresponds to the equivalent cell in the refined list.

Sources
=======

General notes
-------------

Sources can be added to the model using methods of the form
``m.add_*_source()``. For example adding a point source can be done with::

    source = m.add_point_source()

These methods return a handle to the source object, which if captured allow
the user to set and modify the source parameters::

    source = m.add_point_source()
    source.luminosity = lsun
    source.temperature = 10000.
    source.position = (0., 0., 0.)

.. note:: it is also possible to specify the parameters using keyword
          arguments during initalization, e.g.::

              m.add_point_source(luminosity=lsun, temperature=10000.,
                                 position=(0., 0., 0.))

          though this can be longer to read for sources with many arguments.

All sources require a luminosity, given by the ``luminosity`` attribute (or
``luminosity=`` argument), and the emission spectrum can be defined in one of
three ways:

* by specifying a spectrum using the ``spectrum`` attribute (or ``spectrum=``
  argument). The spectrum should either be a tuple of (nu, fnu) or an instance
  of an ``atpy.Table`` with two columns named 'nu' and 'fnu'. For example,
  given a file ``spectrum.txt`` with two columns listing frequency and flux,
  the spectrum can be set using::

    import numpy
    spectrum = np.loadtxt('spectrum.txt', dtype=[('nu', float),
                                                 ('fnu', float)])
    source.spectrum = (spectrum['nu'], spectrum['fnu'])

* by specifying a blackbody temperature using the ``temperature`` attribute
  (or ``temperature=`` argument). This should be a floating point value.

* by using the local dust emissivity if neither a spectrum or temperature are
  specified.

Point Sources
-------------

A point source is defined by a luminosity, a 3-d cartesian position (set to
the origin by default), and a spectrum or temperature. The following examples
demonstrate adding different point sources:

* Set up a 1 solar luminosity 10,000K point source at the origin::

    source = m.add_point_source()
    source.luminosity = lsun  # [ergs/s]
    source.temperature = 10000.  # [K]

* Set up two 0.1 solar luminosity 1,300K point sources at +/- 1 AU in the x
  direction::

    # Set up the first source
    source1 = m.add_point_source()
    source1.luminosity = 0.1 * lsun  # [ergs/s]
    source1.position = (au, 0, 0)  # [cm]
    source1.temperature = 1300.  # [K]

    # Set up the second source
    source2 = m.add_point_source()
    source2.luminosity = 0.1 * lsun  # [ergs/s]
    source2.position = (-au, 0, 0)  # [cm]
    source2.temperature = 1300.  # [K]

* Set up a 10 solar luminosity source at the origin with a spectrum read in
  from a file with two columns giving wav (in microns) and fnu::

    # Use NumPy to read in the spectrum
    import numpy as np
    data = np.loadtxt('spectrum.txt', dtype=[('wav', float), ('fnu', float)])

    # Convert to nu, fnu
    nu = c / (data['wav'] * 1.e-4)
    fnu = data['fnu']

    # Set up the source
    source = m.add_point_source()
    source.luminosity = 10 * lsun  # [ergs/s]
    source.spectrum = (nu, fnu)

Spherical Sources
-----------------

Adding spherical sources is very similar to adding point sources, with the
exception that a radius can be specified::

    source = m.add_spherical_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = rsun  # [cm]
    source.temperature = 10000.  # [K]

It is possible to add limb darkening, using::

    source.limb = True

Spots
-----

Adding spots to a spherical source is straightforward. Spots behave the same
as other sources, requiring a luminosity, spectrum, and additional geometrical
parameters::

    source = m.add_spherical_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = rsun  # [cm]
    source.temperature = 10000.  # [K]

    spot = source.add_spot()
    spot.luminosity = 0.1 * lsun  # [ergs/s]
    spot.longitude = 45.  # [degrees]
    spot.latitude = 30.  # [degrees]
    spot.radius = 5.  # [degrees]
    spot.temperature = 20000.  # [K]

Map Sources
-----------

Map sources are diffuse sources that are defined by a total luminosity, and a
probability distribution map for the emission, defined on the same grid as the
density. For example, if the grid is defined on a 10x10x10 grid, the following
will add a source which emits photons from all cells equally::

    source = m.add_map_source()
    source.luminosity = lsun  # [ergs/s]
    source.map = np.ones((10, 10, 10))

.. note:: The ``map`` array does not need to be normalized.

External sources
----------------

There are two kinds of external illumination sources, spherical and box
sources - the former being more suited to spherical polar grids, and the
latter to cartesian, AMR, and octree grids (there is no cylindrical external
source for cylindrical grids at this time). For example, an external spherical source can be added with::

    source = m.add_external_spherical_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = pc  # [cm]
    source.temperature = 10000.  # [K]

As for point and spherical sources, the position of the center can also be
set, and defaults to the origin. External box sources have a ``bounds`` attribute instead of ``radius`` and ``position``::

    source = m.add_external_box_source()
    source.luminosity = lsun  # [ergs/s]
    source.bounds = [[-pc, pc], [-pc, pc], [-pc, pc]]  # [cm]
    source.temperature = 10000.  # [K]

where the ``bounds`` attribute is given as
``[[xmin, xmax], [ymin, ymax], [zmin, zmax]]``.

Plane parallel sources
----------------------

Finally, it is possible to add circular plane parallel sources (essentially a
circular beam with a given origin and direction)::

    source = m.add_external_spherical_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = rsun  # [cm]
    source.temperature = 10000.  # [K]
    source.position = (au, 0., 0.)  # [cm]
    source.direction = (45., 0.)  # [degrees]

where ``direction`` is a tuple of (theta, phi) that gives the direction of the
beam.

Configuration
=============

To configure the parameters for the model, such as number of photons or number of iterations, the following methods are available::

Number of photons
-----------------

The number of photons to run in various iterations is set using the
following method::

    m.set_n_photons(...)

This method can take the following arguments, which depend on the type of radiation transfer calculations requested:

* ``initial=`` - number of photons per initial iteration to compute the
  specific energy of the dust
* ``imaging=`` - number of photons emitted in the SED/image iteration.
* ``raytracing_sources=`` - number of photons emitted from sources in the
  raytracing iteration
* ``raytracing_dust=`` - number of photons emitted from dust in the raytracing
  iteration
* ``stats=`` - used to determine how often to print out statistics

If computing the radiation transfer in monochromatic mode, the ``imaging`` argument should be replaced by:

* ``imaging_sources=`` - number of photons emitted from sources in the
  SED/image iteration.
* ``imaging_dust=`` - number of photons emitted from dust in the SED/image
  iteration.

.. note:: Only the relevant arguments need to be specified - for example if no
          sources are present, the ``*_sources`` arguments can be ignored,
          while if no dust density grids are present, the ``*_dust`` arguments
          can be ignored.

.. note:: All the required arguments have to be specified in a single call to
          ``set_n_photons``.

Specific Energy calculation
---------------------------

To set the number of initial iterations used to compute the dust specific
energy, use::

    m.set_n_initial_iterations(10)

Raytracing
----------

To enable raytracing, simply use::

    m.set_raytracing(True)

Diffusion
---------

If the model density contains regions of very high density where photons
get trapped or do not enter, one can enable either or both the modified
random walk (MRW; Min et al. 2009, Robitaille et al. 2010) and the partial
diffusion approximation (PDA; Min et al. 2009). The MRW requires a
parameter ``gamma`` which is used to determine when to start using the MRW
(see Min et al. 2009 for more details). By default, this parameter is set
to one. The following examples show how to enable the PDA and MRW respectively:

* Enable the partial diffusion approximation::

    m.set_pda(True)

* Enable the modified random walk, and set the gamma parameter to 2::

    m.set_mrw(True, gamma=2)

Dust sublimation
----------------

To set whether and how to sublimate dust, first the dust file needs to be read in, the sublimation parameters should be set, and the dust object should be passed directly to add_density::

    from hyperion.dust import SphericalDust

    dust = SphericalDust('kmh.hdf5')
    dust.set_sublimation_temperature('fast', temperature=1600)

    m.add_density_grid(density, dust)

The first argument of ``set_sublimation_temperature`` can be ``none`` (dust sublimation does not occur), ``cap`` (temperatures in excess of the one specified will be reset to the one given), ``slow`` (dust with temperatures in excess of the one specified will be gradually destroyed), or ``fast`` (dust with temperatures in excess of the one specified will be immediately destroyed).

Advanced Settings
-----------------

Set the maximum number of photon interactions::

    m.set_max_interactions(100000)

Kill all photons as soon as they are absorbed, in the imaging/SED iteration
(not in the temperature iterations)::

    m.set_kill_on_absorb(True)

Set a minimum temperature to which temperatures below this will be reset::

    m.add_density_grid(density, dust, minimum_temperature=100.)

and in terms of specific energy::

    m.add_density_grid(density, dust, minimum_specific_energy=100.)

Set the number of output bytes per floating point value (4 = 32-bit, 8 = 64-bit)::

    m.set_output_bytes(4)


