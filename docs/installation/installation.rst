============
Installation
============

The easy way
============

The easiest way to install Hyperion and all the dependencies on MacOS X or Linux
is to use ``conda``. If you have never used ``conda`` before, the easiest way to
get started is to install `Miniforge <https://conda-forge.org/miniforge/>`_. You
will need to download the appropriate **Miniforge3-*.sh** file and install it by
running e.g.::

    bash Miniforge3-24.7.1-0-Linux-x86_64.sh

You should replace the name of the installer with the one you have downloaded.

Press 'enter' to view the license, then 'q' to proceed, and type ``yes`` to
accept the license. Then, accept the default location (unless you need to change
this), and finally you likely want to answer ``yes`` to the final question so
that conda gets activated by default. As some versions of Miniforge3 do not
always do this properly, open a new terminal and make sure you explicitly tell
conda to activate the base environment automatically::

     conda config --set auto_activate_base true

Once this is set up, or if you already have a ``conda`` distribution, you can
install Hyperion by opening a new terminal and doing::

    conda install -c conda-forge hyperion

This will install both the Fortran binaries and the Python library for Hyperion
(as well as all the dependencies including MPI and HDF5). That's it! You can
check that the installation works by making sure that the following commands do
not return 'command not found'::

    $ hyperion
    usage: hyperion [-h] [-f] [-m n_cores] input output
    hyperion: error: the following arguments are required: input, output
    $ hyperion_car
    Usage: hyperion_car [-f] input_file output_file

If you see the same as the above, you are all set! If you run into any issues,
please let us know by `opening an issue
<https://github.com/hyperion-rt/hyperion/issues>`_.

The longer way
==============

Fortran dependencies
--------------------

The packages required for the Fortran part of Hyperion are:

* A Fortran compiler. The following compilers are known to work:

  * gfortran 4.3 and later
  * ifort 11 and later
  * pgfortran 11 and above

* `HDF5 <http://www.hdfgroup.org/HDF5/>`_ 1.8.x or 1.10.x with the Fortran bindings

* An MPI installation (e.g. `MPICH2 <http://www.mpich.org/>`_ or `OpenMPI
  <http://www.open-mpi.org/>`_) with the Fortran bindings

Note that in some cases, default installations of HDF5 and MPI packages do not
include support for Fortran - this has to be explicitly enabled.

Due to the variety of operating system versions and Linux distributions, we
can't provide detailed instructions for each one, so if you need help with
installing these dependencies, consider getting help from your local friendly
system administrator.

There are nevertheless a couple of simple cases. First, on Debian-based Linux
distributions (including Ubuntu), you should be able to install these
dependencies with::

    apt-get install libmpich2-dev libhdf5-serial-dev

On Fedora Linux distributions, you should be able to install these
dependencies with::

    yum install hdf5-static mpich-devel libgfortran-static

Python dependencies
-------------------

The packages required for the Python code (in addition to Python itself) are:

* `NumPy <http://www.numpy.org>`_ 1.11 or later
* `Matplotlib <http://matplotlib.org>`_ 1.5 or later
* `h5py <http://www.h5py.org>`_ 2.4 or later
* `Astropy <http://www.astropy.org>`_ 1.2 or later
* `yt <http://yt-project.org/>`_ 3.3 or later

These dependencies will be automatically installed when installing the Python
component of Hyperion if they are not already present.

.. _hyperion_install:

Hyperion
--------

Once you have installed the dependencies, you are ready to install Hyperion!

Download the latest tar file from `here <https://pypi.python.org/pypi/Hyperion/>`_, then expand it with::

    tar xvzf Hyperion-x.x.x.tar.gz
    cd Hyperion-x.x.x

Python module
^^^^^^^^^^^^^

Install the Python module with::

    pip install ".[recommended]"

or::

    pip install ".[recommended]" --user

if you do not have root access. Check that the module installed correctly::

    $ python
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import hyperion
    >>>

and also try typing::

    $ hyperion

in your shell. If you get ``command not found``, you need to ensure that the
scripts installed by Python are in your ``$PATH``. If you do not know where
these are located, check the last line of the install command above, which
should contain something like this::

    changing mode of /Users/tom/Library/Python/2.7/bin/hyperion to 755

The path listed (excluding ``hyperion`` at the end) should be in your
``$PATH``. On Linux systems, this path will often be ``$HOME/.local/bin``.

.. note:: On recent versions of MacOS X, you may encounter the following error
          when trying to install the Python library for Hyperion::

              clang: error: unknown argument: '-mno-fused-madd' [-Wunused-command-line-argument-hard-error-in-future]

          If this is the case, try setting the following environment variables
          before installing it::

              export CFLAGS=-Qunused-arguments
              export CPPFLAGS=-Qunused-arguments


Fortran binaries
^^^^^^^^^^^^^^^^

If you are using HDF5 1.10.x, compile the Fortran code with::

    ./configure
    make
    make install

If you are using HDF5 1.8.x, compile the Fortran code with::

    HYPERION_HDF5_VERSION=18 ./configure
    make
    make install

By default, the binaries will be written to ``/usr/local/bin`` (which will
require you to use ``sudo`` for the last command). If you would prefer to
install to your home directory, you can change this using the ``--prefix``
option to configure, for example::

    ./configure --prefix=$HOME/usr

To check that the Fortran binaries are correctly installed, try typing::

    $ hyperion_sph
    Usage: hyperion input_file output_file

If you get::

    $ hyperion_sph
    hyperion_sph: command not found

then something went wrong in the installation, or the directory to which you
installed the binaries is not in your ``$PATH``. Otherwise, you are all set!
