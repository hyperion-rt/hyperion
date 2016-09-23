============
Installation
============

.. important:: This section contains information on setting up the dependencies for
          Hyperion as well as Hyperion itself. If you have any issues with the
          installation of any of the dependencies or Hyperion, please first
          talk to your system administrator to see if they can help you get
          set up!

The easy way
============

The easiest way to install Hyperion and all the dependencies on MacOS X or Linux
is to use the `Anaconda Python Distribution <https://www.continuum.io/downloads>`_
or `Miniconda <http://conda.pydata.org/miniconda.html>`_. Once you have either
of these set up, you can install Hyperion by simply doing::

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

If you see the same as the above, you are all set!

The longer way
==============

Dependencies
------------

First, you will need to install several dependencies for the Fortran and Python
versions of Hyperion. Choose your own adventure!

.. toctree::
   :maxdepth: 1

   install_debian_ubuntu.rst
   install_fedora.rst
   install_centos_scilinux.rst
   install_linux_nonroot.rst
   install_macosx.rst
   install_full.rst

.. _hyperion_install:

Hyperion
--------

Once you have installed the dependencies as described in one of the sections
above, you are ready to install Hyperion!

Download the latest tar file from `here <https://pypi.python.org/pypi/Hyperion/>`_, then expand it with::

    tar xvzf Hyperion-x.x.x.tar.gz
    cd Hyperion-x.x.x

Python module
-------------

Install the Python module with::

    python setup.py install

or::

    python setup.py install --user

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
----------------

Compile the Fortran code with::

    ./configure
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

CMake build system
^^^^^^^^^^^^^^^^^^

An experimental build system for the Fortran binaries based on CMake is now
available. You can find the detailed instructions on how to use it at the page
:doc:`cmake_build_system`.
