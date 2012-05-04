Using previous geometry and/or quantities
=========================================

In some cases, it can be useful to re-use a geometry and optionally a density
and/or specific energy grid from a previous model in a new model. For example,
one might want to separate the calculation of the temperature/specific energy
from the calculation of the images. Two methods are available for this:
``use_geometry``, and ``use_quantities``.

First, you can re-use the geometry (i.e. the definition of the grid, excluding
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
    m1.write('model1.rtin')

and run the model to produce ``model1.rtout``, then you can create a new model
that makes use of the geometry to set up a model with the same grid, but
different density values and source properties by doing::

    m2 = Model()
    m2.use_geometry('model1.rtout')
    m2.add_density_grid(np.array([[[2.e-10]]]), 'kmh.hdf5')
    s = m2.add_point_source()
    s.luminosity = 0.5 * lsun
    s.temperature = 6000.
    m2.set_n_photons(initial=1000, imaging=1000)
    m2.write('model2.rtin')

Similarly, you can also instruct Hyperion to use the same density grid by
doing::

    m.use_quantities(<filename>, quantities=<quantities to use>,
                     use_minimum_specific_energy=<bool>, use_dust=<bool>)

In this case, the previous model has to be an output file from Hyperion. By
default, both the density and the specific energy from the last iteration of a
previous model are used, as well as the setting for the minimum specific
energy (if set), and the dust properties. For example, in the above example,
you can do::

    m2 = Model()
    m2.use_geometry('model1.rtout')
    m2.use_quantities('model1.rtout')
    s = m2.add_point_source()
    s.luminosity = lsun
    s.temperature = 6000.
    m2.set_n_photons(initial=1000, imaging=1000)
    m2.write('model2.rtin')

to use the density, specific energy, and dust from ``model1.rtout``. Note
however that the specific energy will be recalculated and overwritten unless
you disable the calculation of the specific energy::

    m2 = Model()
    m2.use_geometry('model1.rtout')
    m2.use_quantities('model1.rtout')
    s = m2.add_point_source()
    s.luminosity = lsun
    s.temperature = 6000.
    m2.set_n_initial_iterations(0)  # disable specific energy calculation
    m2.set_n_photons(imaging=1000)  # don't specify initial number of photons
    m2.write('model2.rtin')

Note that you can also use just the density or just the specific energy if you
wish, by using the ``quantities`` argument, e.g.::

    m2.use_quantities('model1.rtout', quantities=['density'])

or::

    m2.use_quantities('model1.rtout', quantities=['specific_energy'])

You can disable using the dust from the previous model (in case you want to
change it)::

    m2.use_quantities('model1.rtout', use_dust=False)

and you can also disable using the minimum specific energy::

    m2.use_quantities('model1.rtout', use_minimum_specific_energy=False)
