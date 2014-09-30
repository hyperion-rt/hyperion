Radiative transfer settings
===========================

.. _`Min et al. 2009`: http://www.aanda.org/index.php?option=com_article&access=bibcode&Itemid=129&bibcode=2009A%2526A...497..155MFUL

.. _`Robitaille 2010`: http://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201015025&Itemid=129

.. _`Robitaille (2011)`: http://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201117150&Itemid=129>

Once the coordinate grid, density structure, dust properties, and luminosity
sources are set up, all that remains is to set the parameters for the
radiation transfer algorithm, including number of photons to use, or whether
to use various optimization schemes.

Number of photons
-----------------

The number of photons to run in various iterations is set using the
following method::

    m.set_n_photons(initial=1000000, imaging=1000000)

where ``initial`` is the number of photons to use in the iterations for
the specific energy (and therefore temperature), and ``imaging`` is the
number of photons for the SED/image calculation, whether using binned
images/SEDs or peeling-off.

In addition, the ``stats=`` argument can be optionally specified to indicate
how often to print out performance statistics (if it is not specified a
sensible default is chosen).

Since the number of photons is crucial to produce good quality results, you
can read up more about setting sensible values at :doc:`photon_numbers`.

.. _convergence:

Specific energy calculation
---------------------------

To set the number of initial iterations used to compute the dust specific
energy, use e.g.::

    m.set_n_initial_iterations(5)

Note that this can also be zero, in which case the temperature is not solved,
and the radiative transfer calculation proceeds to the image/SED calculation
(this is useful for example if one is making images at wavelengths where
thermal emission is negligible, or if a specific energy/temperature was
specified as input).

It is also possible to tell the radiative transfer algorithm to exit these
iterations early if the specific energy has converged. To do this, use::

    m.set_convergence(True, percentile=100., absolute=0., relative=0.)

where the boolean value indicates whether to use convergence detection
(``False`` by default), and ``percentile``, ``absolute``, and ``relative``
arguments are explained in more detail in Section 2.4 of `Robitaille (2011)`_.
For the benchmark problems of that paper, the values were set to::

    m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

which are reasonable starting values. Note that if you want to use convergence
detection, you should make sure that the value for
``set_n_initial_iterations`` is not too small, otherwise the calculation might
stop before converging. When running the main Hyperion code, convergence
statistics are printed out, and it is made clear when the specific energy has
converged.

.. _initial_specific_energy:

Initial and additional specific energy
--------------------------------------

Another option that is related to the specific energy is
:meth:`~hyperion.model.Model.set_specific_energy_type`. This is used to control
how any specific energy passed to
:meth:`~hyperion.model.Model.add_density_grid` is used. By default, the
specific energy specified is the *initial* specific energy used, and if the
number of temperature iterations is not zero (see :ref:`convergence`) this
specific energy gets replaced with the self-consistently calculated one in
later iterations. If instead you want this specific energy to be *added* to the
self-consistently computed one after each iteration, you can set::

    m.set_specific_energy_type('additional')

This can be used for example if you need to take into account an additional
source of heating that cannot be modelled by Hyperion.

Raytracing
----------

To enable raytracing (for source and dust emission, but not scattering),
simply use::

    m.set_raytracing(True)

This algorithm is described in Section 2.6.3 of `Robitaille (2011)`_. If raytracing is used, you will need to add the ``raytracing_sources`` and ``raytracing_dust`` arguments to the call to ``set_n_photons``, i.e.::

    m.set_n_photons(initial=1000000, imaging=1000000,
                    raytracing_sources=1000000, raytracing_dust=1000000)

.. _diffusion:

Diffusion
---------

If the model density contains regions of very high density where photons get
trapped or do not enter, one can enable the modified random walk (MRW; `Min et
al. 2009`_, `Robitaille 2010`_) in order to group many photon interactions
into one. The MRW requires a parameter ``gamma`` which is used to determine
when to start using the MRW (see `Min et al. 2009`_ for more details). By
default, this parameter is set to ``1``. The following example shows how to
enable the modified random walk, and set the gamma parameter to ``2``::

    m.set_mrw(True, gamma=2.)

In some cases (such as protoplanetary disks) very optically thick regions do
not receive any radiation. In cases where the temperature in these regions is
important, one can use the partial diffusion approximation (PDA; `Min et al.
2009`_) to solve the diffusion equation over the grid and find the missing
temperatures::

    m.set_pda(True)

Note however that if more than 10,000 cells have low photon counts and require
the PDA, this can be **very** slow, so this option is only recommended in
cases where you know it is absolutely needed. In most cases, if photons cannot
reach inside certain cells, these cells are unlikely to be contributing a
significant amount of flux to SEDs or images.

Dust sublimation
----------------

To set whether and how to sublimate dust, first the dust file needs to be read
in (or initialized in the script), the sublimation parameters should be set,
and the dust object should be passed directly to add_density::

    from hyperion.dust import SphericalDust

    d = SphericalDust('kmh.hdf5')
    d.set_sublimation_temperature('fast', temperature=1600.)

    m.add_density_grid(density, d)

The first argument of ``set_sublimation_temperature`` can be ``none`` (dust
sublimation does not occur), ``cap`` (temperatures in excess of the one
specified will be reset to the one given), ``slow`` (dust with temperatures in
excess of the one specified will be gradually destroyed), or ``fast`` (dust
with temperatures in excess of the one specified will be immediately
destroyed). For more information, see Section 2.7.3 of `Robitaille (2011)`_.

.. _sample_sources_evenly:

Multiple sources
----------------

By default, the number of photons emitted is proportional to the luminosity
of the sources, so in cases where several sources with very different
luminosities are included in the models, some sources might be
under-sampled. In some cases, this will not be a problem, but in some cases
you may want to emit equal numbers of photons from each source instead. For
example, if you have two sources that have a bolometric luminosity that is
different by a factor of 100, but at the specific wavelength you are
interested in they have the same flux, then you will probably want equal
numbers of photons for both sources. You can enable this with::

    m.set_sample_sources_evenly(True)

Note that this does not yet cause point sources within
:ref:`point-source-collections` to be evenly sampled. For the purposes of this
option, a point source collection counts as a single source.

Outputting physical quantities
------------------------------

It is possible to write out a number of physical arrays for each iteration, or
just the last iteration. To do this, you will need to set the parameters in
``Models.conf.output``::

    # Density
    m.conf.output.output_density = 'last'

    # Density difference (shows where dust was destroyed)
    m.conf.output.output_density_diff = 'none'

    # Energy absorbed (using pathlengths)
    m.conf.output.output_specific_energy = 'last'

    # Number of unique photons that passed through the cell
    m.conf.output.output_n_photons = 'last'

Each value can be set to ``all`` (output all iterations), ``last`` (output
only after last iteration), or ``none`` (do not output). The default is to
output only the last iteration of ``specific_energy``. To find out how to view
these values, see :doc:`../postprocessing/postprocessing`

Advanced parameters
-------------------

There are a number of more advanced parameters to control the radiative
transfer, but since they are not essential initially, they are described in
the :doc:`../advanced/advanced_conf` section.
