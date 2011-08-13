==============
Running models
==============

Once an .rtin file has been created (see :ref:`setup`), the model can be run using the compiled Fortran code. Note that the model can be run on a different computer/cluster to the computer on which is was set up, because the ``.rtin`` files are portable.

The easiest way to run a model is to invoke the ``hyperion`` command-line utility, specifying the input and output file (we use the ``.rtout`` extensions for output files)::

    hyperion model.rtin model.rtout

``hyperion`` is in fact a wrapper to grid specific binaries:

* ``hyperion_car`` for cartesian grids
* ``hyperion_cyl`` for cylindrical polar grids
* ``hyperion_sph`` for spherical polar grids
* ``hyperion_amr`` for AMR grids
* ``hyperion_oct`` for Oct-tree grids

These binaries can be called directly if necessary.

To use the parallel version of the code, use the relevant binary, with the ``_mpi`` suffix appended, and launch it using the command relevant to your MPI installation, for example::

    mpirun -n 128 hyperion_car_mpi model.rtin model.rtout

This can also be ``mpiexec`` or ``openmpirun`` or ``openmpiexec`` depending on your MPI installation. Note that there is no wrapper to select the correct grid for the MPI-enabled code.

It is also possible to run the serial version of the code directly from the set-up script, by doing::

    ...
    m.write()
    m.run()

