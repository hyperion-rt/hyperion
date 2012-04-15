Radiative transfer settings
===========================

To configure the parameters for the model, such as number of photons or number
of iterations, the following methods are available::

Number of photons
-----------------

The number of photons to run in various iterations is set using the
following method::

    m.set_n_photons(...)

This method can take the following arguments, which depend on the type of
radiation transfer calculations requested:

* ``initial=`` - number of photons per initial iteration to compute the
  specific energy of the dust
* ``imaging=`` - number of photons emitted in the SED/image iteration.
* ``raytracing_sources=`` - number of photons emitted from sources in the
  raytracing iteration
* ``raytracing_dust=`` - number of photons emitted from dust in the raytracing
  iteration
* ``stats=`` - used to determine how often to print out statistics

If computing the radiation transfer in monochromatic mode, the ``imaging``
argument should be replaced by:

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
(see Min et al. 2009 for more details). By default, this parameter is set to
one. The following examples show how to enable the PDA and MRW respectively:

* Enable the partial diffusion approximation::

    m.set_pda(True)

* Enable the modified random walk, and set the gamma parameter to 2::

    m.set_mrw(True, gamma=2)

Dust sublimation
----------------

To set whether and how to sublimate dust, first the dust file needs to be read
in, the sublimation parameters should be set, and the dust object should be
passed directly to add_density::

    from hyperion.dust import SphericalDust

    dust = SphericalDust('kmh.hdf5')
    dust.set_sublimation_temperature('fast', temperature=1600)

    m.add_density_grid(density, dust)

The first argument of ``set_sublimation_temperature`` can be ``none`` (dust
sublimation does not occur), ``cap`` (temperatures in excess of the one
specified will be reset to the one given), ``slow`` (dust with temperatures in
excess of the one specified will be gradually destroyed), or ``fast`` (dust
with temperatures in excess of the one specified will be immediately
destroyed).

Physical quantity output
------------------------

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
output only the last iteration of ``specific_energy``.

Advanced Settings
-----------------

Set the maximum number of photon interactions::

    m.set_max_interactions(100000)

Kill all photons as soon as they are absorbed, in the imaging/SED iteration
(not in the temperature iterations)::

    m.set_kill_on_absorb(True)

Set the number of output bytes per floating point value (4 = 32-bit, 8 =
64-bit)::

    m.set_output_bytes(4)

For the :doc:`model` class, one can set a minimum temperature to which
temperatures below this will be reset::

    m.add_density_grid(density, dust, minimum_temperature=100.)

and in terms of specific energy::

    m.add_density_grid(density, dust, minimum_specific_energy=100.)

For the :doc:`analytical_yso_model`, this can be done with::

    m.set_minimum_temperature(100.)
