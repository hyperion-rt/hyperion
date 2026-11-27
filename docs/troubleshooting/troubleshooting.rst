.. _report: http://www.github.com/hyperion-rt/hyperion/issues

Troubleshooting
===============

This page describes problems that have been reported when installing or
running Hyperion, and how to solve them. If you run into a problem that is not
listed here, please `report`_ it.

A run stops with "output file appears to be corrupt"
----------------------------------------------------

When the Fortran code stops with an error, the ``hyperion`` command prints::

    Run did not complete successfully: output file appears to be corrupt

and ``Model.run`` raises::

    SystemExit: An error occurred, and the run did not complete

These messages only mean that the run did not finish: the actual error is
printed earlier in the output (or in the log file, if one was given to
``Model.run``), so look further up for the first error message.

"cannot read files made with the Python module before version 0.8.7"
--------------------------------------------------------------------

If the Fortran code stops with::

    ERROR: cannot read files made with the Python module before version 0.8.7

(or the equivalent message for dust files) even though you are using a recent
version of the Python module, check the version reported by::

    python -c "import hyperion; print(hyperion.__version__)"

If this gives a version such as ``0.1.dev1628``, Hyperion was installed from a
git repository that does not include the version tags, which are used to
determine the version. This typically happens when installing from a fork. To
fix this, fetch the tags from the main repository and re-install the Python
module::

    git fetch --tags https://github.com/hyperion-rt/hyperion.git
    pip install .

MPI fails to start with errors mentioning UCX
---------------------------------------------

When running in parallel, MPI may fail to start with errors such as::

    UCX  ERROR ibv_create_srq() failed: Operation not supported
    Abort(...) on node 0: Fatal error in internal_Init: Other MPI error, error stack:
    MPI_Init(argc=(nil), argv=(nil)) failed
    ...
    init_worker(38)...............:  ucx function returned with failed status

This is not specific to Hyperion: it happens with MPI libraries that use the
UCX communication library (such as the MPICH package in Ubuntu 24.04) on
machines with network devices that appear to support InfiniBand but do not
support all its features, which is common on cloud virtual machines. When all
the processes run on a single machine, you can restrict UCX to shared memory
and TCP by setting::

    export UCX_TLS=self,sm,tcp

before running Hyperion. On a cluster with a working InfiniBand network, do
not set this variable, since it would prevent MPI from using the InfiniBand
network between nodes.

Open MPI reports that there are not enough slots available
----------------------------------------------------------

When running with Open MPI, using more processes than there are cores on the
machine fails with::

    There are not enough slots available in the system to satisfy the 8
    slots that were requested by the application

since Open MPI does not run more than one process per core by default. Either
use fewer processes, or allow Open MPI to run more processes than there are
cores by setting the following before running Hyperion::

    export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe

for Open MPI 5, or::

    export OMPI_MCA_rmaps_base_oversubscribe=true

for Open MPI 4. Note that running more processes than there are cores will not
make the code run faster.

Firewall prompts on Macs when running with MPI
----------------------------------------------

When running the parallel code on a Mac with the firewall enabled, you will
typically be asked once::

    Do you want the application "prte" to accept incoming network connections?

``prte`` is the program that Open MPI (version 5 and later) uses to launch the
processes, and it listens for connections from them. It is safe to allow this,
and macOS remembers the answer until Open MPI is updated.

In addition, each of the ``hyperion_*_mpi`` processes opens network ports on
all network interfaces to communicate with the other processes, which can
result in one prompt per process on every run, since binaries compiled locally
are not recognized as the same application once they are recompiled. When all
the processes run on the same machine, they can instead communicate through
shared memory, in which case they do not open any network ports. With Open
MPI, this can be done by setting the following in ``~/.hyperionrc``::

    [mpi]
    command = mpirun --mca btl self,sm

(see :doc:`../running/running`). The same options can be passed to ``mpirun``
when calling the Fortran binaries directly. You can check which programs are
listening for network connections during a run with::

    lsof -nP -iTCP -sTCP:LISTEN | grep -i -E "prte|hyperion"
