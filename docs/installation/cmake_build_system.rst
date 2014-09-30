==============================
Alternative CMake build system
==============================

An experimental build system for the Fortran binaries based on
`CMake <http://www.cmake.org/>`_ is now available. The two key
advantages with respect to the default build system are:

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

.. NOTE::
   Normally in order to change the compiler it will be necessary
   to completely erase the contents of the build directory and start from
   scratch. This is not necessary when changing other CMake variables
   such as those discussed below.
   The compiler variable is special because CMake uses it as
   a starting point to detect and setup the compilation environment.

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
