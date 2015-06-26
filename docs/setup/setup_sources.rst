Luminosity sources
==================

General notes
-------------

Sources can be added to the model using methods of the form
``m.add_*_source()``. For example adding a point source can be done with::

    source = m.add_point_source()

These methods return a source 'object' that can be used to set and modify the
source parameters::

    source = m.add_point_source()
    source.luminosity = lsun
    source.temperature = 10000.
    source.position = (0., 0., 0.)

.. note:: It is also possible to specify the parameters using keyword
          arguments during initialization, e.g.::

              m.add_point_source(luminosity=lsun, temperature=10000.,
                                 position=(0., 0., 0.))

          though this can be longer to read for sources with many arguments.

All sources require a luminosity, given by the ``luminosity`` attribute (or
``luminosity=`` argument), and the emission spectrum can be defined in one of
three ways:

* by specifying a spectrum using the ``spectrum`` attribute (or ``spectrum=``
  argument). The spectrum should either be a ``(nu, fnu)`` pair or an instance
  of an ``atpy.Table`` with two columns named ``'nu'`` and ``'fnu'``. For
  example, given a file ``spectrum.txt`` with two columns listing frequency
  and flux, the spectrum can be set using::

    import numpy
    spectrum = np.loadtxt('spectrum.txt', dtype=[('nu', float),
                                                 ('fnu', float)])
    source.spectrum = (spectrum['nu'], spectrum['fnu'])

* by specifying a blackbody temperature using the ``temperature`` attribute
  (or ``temperature=`` argument). This should be a floating point value.

* by using the local dust emissivity if neither a spectrum or temperature are
  specified.

.. note:: By default, the number of photons emitted is proportional to the
          luminosity, so in cases where several sources with very different
          luminosities are included in the models, some sources might be
          under-sampled. You can instead change the configuration to emit
          equal number of photons from all sources -
          see :ref:`sample_sources_evenly` for more details.

Point sources
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
  from a file with two columns giving wavelength (in microns) and
  monochromatic flux::

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

.. note:: Regardless of the grid type, the coordinates for the sources should
          always be specified in cartesian coordinates, and in the order
          ``(x, y, z)``.

If you want to set up many point sources (for example for a galaxy model) you
may instead want to consider using a `Point source collections`_.

Spherical sources
-----------------

Adding spherical sources is very similar to adding point sources, with the
exception that a radius can be specified::

    source = m.add_spherical_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = rsun  # [cm]
    source.temperature = 10000.  # [K]

It is possible to add limb darkening, using::

    source.limb = True

Spots on spherical sources
--------------------------

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

Diffuse sources
---------------

Diffuse sources are defined by a total luminosity, and a probability
distribution map for the emission, defined on the same grid as the density.
For example, if the grid is defined on a 10x10x10 grid, the following will add
a source which emits photons from all cells equally::

    source = m.add_map_source()
    source.luminosity = lsun  # [ergs/s]
    source.map = np.ones((10, 10, 10))

By default, if no spectrum or temperature is provided, photons will be emitted
using the local emissivity of the dust. However, you can also specify either a
temperature or a spectrum as for `Point Sources`_ and `Spherical Sources`_, e.g::

    source.temperature = 10000.  # [K]

or::

    source.spectrum = (nu, fnu)

.. note:: The ``map`` array does not need to be normalized.

External sources
----------------

There are two kinds of external illumination sources, spherical and box
sources - the former being more suited to spherical polar grids, and the
latter to cartesian, AMR, and octree grids (there is no cylindrical external
source for cylindrical grids at this time). In both cases, photons are emitted
inwards isotropically. For example, an external spherical source can be added
with::

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

See :doc:`../tutorials/howto_scaling_isrf` for information on setting the
luminosity correctly in order to reproduce a given intensity field.

.. note:: Even though these sources are referred to as 'external', they have
          to be placed inside the outermost walls of the grid. The sources are
          not box-shared source or spherical source that can be placed outside
          the grid, but rather sources that emit inwards instead of outwards,
          making it possible to simulate an external radiation field.

Plane parallel sources
----------------------

Finally, it is possible to add circular plane parallel sources (essentially a
circular beam with a given origin and direction)::

    source = m.add_plane_parallel_source()
    source.luminosity = lsun  # [ergs/s]
    source.radius = rsun  # [cm]
    source.temperature = 10000.  # [K]
    source.position = (au, 0., 0.)  # [cm]
    source.direction = (45., 0.)  # [degrees]

where ``direction`` is a tuple of (theta, phi) that gives the direction of the
beam.

.. _point-source-collections:

Point source collections
------------------------

In cases where you want to set up more than a few dozen point sources, it may
be worth instead using a point source collection, which can contain an
arbitrary number of point sources with different luminosities, and a common
temperature or spectrum. To add a point source collection, use e.g.::

    source = m.add_point_source_collection()

The attributes are the same as for the `Point Sources`_ but the
``source.luminosity`` attribute should be set to an array with as many elements
as sources, and the ``source.position`` attribute should be set to a 2-d array
where the first dimension matches ``source.luminosity``, and with 3 elements in
the second dimension (x, y, and z). The following example shows how to set up
1000 random point sources with random positions from -1au to 1au in all
directions, and with random luminosities between 0 and lsun::

    N = 1000
    x = np.random.uniform(-1., 1, N) * au
    y = np.random.uniform(-1., 1, N) * au
    z = np.random.uniform(-1., 1, N) * au

    source = m.add_point_source_collection()
    source.luminosity = np.random.random(N) * lsun
    source.position = np.vstack([x, y, z]).transpose()
    source.temperature = 6000.

In terms of photon sampling, a point source collection acts as a single source
with a luminosity given by the sum of the components - so if you have one point
source collection and one spherical source with the same total luminosity, the
number of photons will be evenly split between the two. Within the point source
collection, the number of photons is split according to luminosity.
