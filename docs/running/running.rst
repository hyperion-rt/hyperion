.. _report: http://www.github.com/hyperion-rt/hyperion/issues

Running models
==============

Using the `hyperion` command-line wrapper
-----------------------------------------

Once an .rtin file has been created (see :doc:`../setup/setup`), the model can
be run using the compiled Fortran code. Note that the model can be run on a
different computer/cluster to the computer on which is was set up, because the
``.rtin`` files are portable.

The easiest way to run a model is to invoke the ``hyperion`` command-line
utility, specifying the input and output file (we use the ``.rtout``
extensions for output files)::

    hyperion model.rtin model.rtout

.. note:: ``hyperion`` is a command-line Python wrapper around the Fortran
          binaries that gets installed with the Python Hyperion library.

To run Hyperion in parallel, you can use::

    hyperion -m <n_processes> model.rtin model.rtout

where ``<n_processes>`` is the number of processes to run in parallel (does not need to equal the number of cores in the computer or cluster). For example, to run the code over 24 processes, you can use::

    hyperion -m 24 model.rtin model.rtout

This may not work with all MPI installations. If you have issues, see :doc:`../troubleshooting/troubleshooting` for known problems and the next section on calling the Fortran binaries directly (and `report`_ the issue).

The ``hyperion`` command normally launches the parallel code with ``mpirun``
or ``mpiexec``, whichever is found first. If you need to customize the
launcher (for example to pass additional options to it), you can set the
``command`` option in the ``[mpi]`` section of the ``~/.hyperionrc`` file -
the number of processes and binary name are appended to it::

    [mpi]
    command = mpirun

MPI popups and firewall prompts on Mac
--------------------------------------

When running the parallel code on a Mac, you may get one popup per process
asking whether the ``hyperion_*_mpi`` application should be allowed to accept
incoming network connections. This is because each MPI process opens network
sockets on all network interfaces, and because locally compiled binaries are
not signed with a stable identity, so the macOS firewall cannot remember your
answer between runs (or between rebuilds of the binaries).

Since all the communication is local when running on a single machine, the
cleanest solution is to keep MPI off the external network interfaces
entirely, in which case no prompts appear at all. With OpenMPI, this can be
done by setting the following in ``~/.hyperionrc``::

    [mpi]
    command = mpirun --mca btl self,sm --mca oob_tcp_if_include lo0 --mca btl_tcp_if_include lo0

and with MPICH::

    [mpi]
    command = mpiexec -iface lo0

The same options can be passed directly to ``mpirun``/``mpiexec`` when
calling the Fortran binaries directly as described in the next section.

Calling the Fortran binaries directly
-------------------------------------

``hyperion`` is in fact a wrapper to grid specific binaries:

* ``hyperion_car`` and ``hyperion_car_mpi`` for cartesian grids
* ``hyperion_cyl`` and ``hyperion_cyl_mpi`` for cylindrical polar grids
* ``hyperion_sph`` and ``hyperion_sph_mpi`` for spherical polar grids
* ``hyperion_amr`` and ``hyperion_amr_mpi`` for AMR grids
* ``hyperion_oct`` and ``hyperion_oct_mpi`` for Oct-tree grids

These binaries can be called directly instead of the ``hyperion`` wrapper. For example, to run a model with a cartesian grid in serial, you would use::

    hyperion_car model.rtin model.rtout

To use the parallel version of the code, use the relevant binary, with the ``_mpi`` suffix appended, and launch it using the command relevant to your MPI installation, for example::

    mpirun -n 128 hyperion_car_mpi model.rtin model.rtout

This can also be ``mpiexec`` or ``openmpirun`` or ``openmpiexec`` depending on your MPI installation.

Running the model from the Python scripts
-----------------------------------------

It is also possible to run the serial version of the code directly from the set-up script, by doing::

    m = Model()
    ...
    m.write('model.rtin')
    m.run('model.rtout')

To run in parallel, simply do::

    m.run('model.rtout', mpi=True, n_processes=<n_processes>)

As for the ``hyperion`` command-line wrapper, this may not work with all MPI installations.

Overwriting existing output
---------------------------

By default, if the output file already exists, a confirmation message is shown::

    WARNING: File exists:  test.rtout
    The following command will be run:  rm test.rtout
    Do you wish to continue? (y/n)

However, this is not always desirable (for example when submitting jobs to
clusters). To overwrite an existing output file, then use the ``-f`` option
when calling ``hyperion`` or on of the ``hyperion_*`` commands::

    $ hyperion -f input output

of use the ``overwrite=True`` argument when using ``Model.run``::

    m.run('model.rtout', overwrite=True)

Note that ``Model.run`` does not ask for confirmation, and instead raises an
error if the output file already exists and ``overwrite=True`` is not set.

Debugging output
----------------

To turn on all the internal debugging output of the Fortran code, use the
``-d`` option when calling ``hyperion`` or one of the ``hyperion_*``
commands::

    $ hyperion -d input output

Note that this prints messages every time a photon crosses a cell wall, so
the output is very large and is only useful for models with a very small
number of photons. When running in parallel, the output from the different
processes will be interleaved.
