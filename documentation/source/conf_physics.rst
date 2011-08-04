Physical Arrays
===============

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