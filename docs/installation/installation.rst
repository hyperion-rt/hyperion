============
Installation
============

.. important:: This section contains information on setting up the dependencies for
          Hyperion as well as Hyperion itself. If you have any issues with the
          installation of any of the dependencies or Hyperion, please first
          talk to your system administrator to see if they can help you get
          set up!

Dependencies
============

First, you will need to install several dependencies. Please follow the
instructions at the following pages to make sure that you have all the
dependencies installed.

.. toctree::
   :maxdepth: 1

   fortran_dependencies.rst
   python_dependencies.rst

.. note:: For instructions for specific computer clusters, see the :ref:`specific` instead, then proceed to the instructions for installing Hyperion below.

.. _hyperion_install:

Hyperion
========

Download the latest tar file from `here <https://pypi.python.org/pypi/Hyperion/>`_, then expand it with::

    tar xvzf hyperion-x.x.x.tar.gz
    cd hyperion-x.x.x

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
``$PATH``.

Fortran binaries
----------------

Compile the Fortran code with::

    ./configure
    make
    make install

By default, the binaries will be written to ``/usr/local/bin`` (which will
require you to use ``sudo`` for the last command), but you can change this
using the ``--prefix`` option to configure, for example::

    ./configure --prefix=/usr/local/hyperion

or::

    ./configure --prefix=$HOME/usr

To check that the Fortran binaries are correctly installed, try typing::

    $ hyperion_sph
    Usage: hyperion input_file output_file

If you get::

    $ hyperion_sph
    hyperion_sph: command not found

then something went wrong in the installation, or the directory to which you
installed the binaries is not in your ``$PATH``. Otherwise, you are all set!

Alternative CMake build system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An experimental build system based on `CMake <http://www.cmake.org/>`_ is
now available. The two key advantages with respect to the default
build system are:

* support for parallel builds,
* support for dependency tracking (i.e., if a single file in the Hyperion
  Fortran source is changed, only that file is recompiled).

At the present time, the CMake build system has been tested only on Linux
and OSX with the GNU and intel compilers. Support for other platforms
will be added in the future.

After having installed CMake, as a first step we are going to create a
``build_cmake`` directory in the Hyperion root directory::

    mkdir build_cmake

We are going to perform an out-of-source build: all the files generated
during the compilation of Hyperion will be kept inside the ``build_cmake``
directory (so that it easy to start from scratch by simply erasing
the build directory).

The next step is the invocation of CMake from the build dir::

    cd build_cmake
    cmake ../

CMake will try to identify the Fortran compiler according to some
platform-dependent heuristics. If the detected compiler is not the
desired one, it is possible to set a custom compiler via the ``FC``
environment variable. For instance, from bash::

    FC=ifort cmake ../

or::

    FC=h5pfc cmake ../

<<<<<<< HEAD
.. NOTE::
   Normally in order to change the compiler it will be necessary
   to completely erase the contents of the build directory and start from
   scratch. This is not necessary when changing other CMake variables
   such as those discussed below.
   The compiler variable is special because CMake uses it as
   a starting point to detect and setup the compilation environment.
=======
Note that normally in order to change the compiler it will be necessary
to completely erase the contents of the build directory and start from
scratch.
>>>>>>> 18f2bbf... Support for ifort and documentation.

CMake will try to locate Hyperion's dependencies (HDF5, MPI)
automatically. This usually works fine on Linux systems (where the
installation paths are more or less standardised, especially when
relying on the system's package manager), but on OSX systems it might be
necessary to point CMake to the correct paths. This can be easily done
via the text-based CMake GUI, called ``ccmake``::

    ccmake ../

After pressing the letter ``t`` to enter advanced mode, it will be possible
to set variables such as ``HDF5_hdf5_hl_LIBRARY_RELEASE`` to the correct
paths on your system.

You might notice that there are other interesting options selectable from
the ``ccmake`` GUI. For instance, the variable ``CMAKE_BUILD_TYPE`` selects
the type of build to perform. The default build type is ``Release``; while
developing, the ``Debug`` mode could be more useful. Note that any
variable visualised in the GUI can also be set from the command line
(see below for some examples). Once all the options have been set, we can run
CMake again from the GUI by pressing ``c`` twice, followed by ``g``.

The output of a successful CMake run is a set of Makefiles that can now be
used via the standard ``make`` command (always from the build directory)::

    make
    make install

It is possible to run the compilation in parallel via the ``-j`` switch, e.g.,::

     make -j 8

Note that if you edit an existing Fortran file in the Hyperion source tree,
you do not need to re-run cmake. Invoking ``make`` as usual will be enough.

Complete CMake command-line examples
""""""""""""""""""""""""""""""""""""

Minimal default configuration::

    cmake ../

Override the compiler::

    FC=h5pfc cmake ../

Override the compiler and set ``Debug`` mode::

    FC=h5pfc cmake -DCMAKE_BUILD_TYPE=Debug ../

Override the compiler, set ``Debug`` mode and set custom installation prefix::

    FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/usr/local ../
