Re-using previous models
========================

In some cases, it can be useful to re-use the geometry, physical quantities,
sources, or other parameters, from a previous model. For example, one might
want to separate the calculation of the temperature/specific energy from the
calculation of the images. A number of methods are available for this, and are
described in the sections below.

Re-using a whole model
----------------------

The simplest case is that you may want to read in a previous model, modify it,
and write it out/run it again. The easiest way to do this is to use the
:meth:`~hyperion.model.Model.read` method::

    m = Model.read('some_model.rtin')

Once the model has been read in, it is possible to modify any of the
parameters, add more density grids, sources, and change the parameters.

It is also possible to read in a model from an output file. In this case, what
is read in are the initial parameters/settings/quantities for the model. If
you would like to use the final specific energy (and optionally density if
available), you can call ``read`` with the ``only_initial=`` argument set to
``False``:

    m = Model.read('some_model.rtout', only_initial=False)

If, instead of reading in the whole model, you want to re-use only certain
aspects of a previous model, see the following sections.

Geometry
--------

It is possible to re-use the geometry (i.e. the definition of the grid, excluding
quantities such as density of specific energy) from either a previous input or
output file by doing::

    m = Model()
    m.use_geometry(<filename>)

For example, if you create a model with::

    m1 = Model()
    m1.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m1.add_density_grid(np.array([[[1.e-10]]]), 'kmh.hdf5')
    s = m1.add_point_source()
    s.luminosity = lsun
    s.temperature = 6000.
    m1.set_n_photons(initial=1000, imaging=1000)
    m1.write('model_1.rtin')

and run the model to produce ``model_1.rtout``, then you can create a new model
that makes use of the geometry to set up a model with the same grid, but
different density values and source properties by doing::

    m2 = Model()
    m2.use_geometry('model_1.rtout')
    m2.add_density_grid(np.array([[[2.e-10]]]), 'kmh.hdf5')
    s = m2.add_point_source()
    s.luminosity = 0.5 * lsun
    s.temperature = 6000.
    m2.set_n_photons(initial=1000, imaging=1000)
    m2.write('model_2.rtin')

The :meth:`~hyperion.model.Model.use_geometry` method can take either a
previous input or output file. See :meth:`~hyperion.model.Model.use_geometry`
for more information.

Quantities
----------

Similarly, you can also instruct Hyperion to use the same density grid by
doing::

    m.use_quantities(<filename>, quantities=<quantities to use>)

As for the geometry, the file can be a previous input or output file from
Hyperion. If an input file, then by default the previous input density (and
optionally specific energy) will be used, whereas if an output file, then by
default the final specific energy and the initial density will be used.
By default, this will also read in the minimum specific energy requested for
the grids, and the dust properties.

For example, in the example mentioned in `Geometry`_ you can do::

    m2 = Model()
    m2.use_geometry('model_1.rtout')
    m2.use_quantities('model_1.rtout')
    s = m2.add_point_source()
    s.luminosity = lsun
    s.temperature = 6000.
    m2.set_n_photons(initial=1000, imaging=1000)
    m2.write('model_2.rtin')

to use the density, specific energy, and dust from ``model_1.rtout``. If you
want to keep the specific energy as-is and avoid computing it further, you
should make sure that you disable the calculation of the specific energy::

    m2 = Model()
    m2.use_geometry('model_1.rtout')
    m2.use_quantities('model_1.rtout')
    s = m2.add_point_source()
    s.luminosity = lsun
    s.temperature = 6000.
    m2.set_n_initial_iterations(0)  # disable specific energy calculation
    m2.set_n_photons(imaging=1000)  # don't specify initial number of photons
    m2.write('model_2.rtin')

Note that you can also use just the density or just the specific energy if you
wish, by using the ``quantities`` argument, e.g.::

    m2.use_quantities('model_1.rtout', quantities=['density'])

or::

    m2.use_quantities('model_1.rtout', quantities=['specific_energy'])

In the case where quantities are being read from an output file, you can also
explicitly request that only the *input* quantities be read in::

    m2.use_quantities('model_1.rtout', only_initial=True)

You can disable using the dust from the previous model (in case you want to
change it)::

    m2.use_quantities('model_1.rtout', use_dust=False)

and you can also disable using the minimum specific energy::

    m2.use_quantities('model_1.rtout', use_minimum_specific_energy=False)

If you are computing a model where the density is changing from one iteration
to the next (for example due to dust sublimation), and if you want to use the
final density, you will need to make sure that you run the initial model with
the option to output the density at the last iteration::

    m1.conf.output.output_density = 'last'

Finally, by default the behavior of
:meth:`~hyperion.model.Model.use_quantities` is to read in the data, so that
it can be modified, but if you do not plan to modify the density, specific
energy, or dust properties, you can also simply link to the previous quantities by doing::

    m2.use_quantities(..., copy=False)

For more information, see :meth:`~hyperion.model.Model.use_quantities`.

Sources
-------

You can import sources from a previous input or output file with::

    m2.use_sources('model_1.rtout')

This will read in the sources, and you can then modify them if needed, or add
new ones to the model. For more information, see :meth:`~hyperion.model.Model.use_sources`.

Configuration
-------------

Several methods are available to read in the image/SED configuration, runtime
parameters, and output parameters from a previous model::

    m1.use_image_config(filename)
    m1.use_run_config(filename)
    m1.use_output_config(filename)

As for the `Sources`_, it is then possible to modify these parameters,
and optionally add new images. For more information, see
:meth:`~hyperion.model.Model.use_image_config`,
:meth:`~hyperion.model.Model.use_run_config`, and
:meth:`~hyperion.model.Model.use_output_config`.
