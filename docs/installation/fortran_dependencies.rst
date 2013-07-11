.. _MPI: http://www.mpich.org/downloads/
.. _HDF5 downloads: http://www.hdfgroup.org/ftp/HDF5/current/src/


.. _fortrandep:

=========================
Fortran code dependencies
=========================

Summary of dependencies
=======================

The packages required for the Fortran code are:

* A recent Fortran compiler. The following compilers are known to work:

  * gfortran 4.3 and later
  * ifort 11 and later
  * pgfortran 11 and above

* `HDF5 <http://www.hdfgroup.org/HDF5/>`_ 1.8.5 or later with the Fortran bindings

* An MPI installation (e.g. `MPICH2 <http://www.mpich.org/>`_ or `OpenMPI <http://www.open-mpi.org/>`_) with the Fortran bindings

Note that often, default installations of HDF5 and MPI packages do not include support for Fortran - this has to be explicitly enabled as described below.

Fortran compiler
================

The first dependency is a Fortran compiler. In addition to commercial
compilers (e.g. ``ifort``, ``pgfortran``, ...), there are a couple of free
ones, the most common of which is ``gfortran``. If you don't already have a
compiler installed, you can install ``gfortran`` via package managers on Linux
machines, or from MacPorts or binary installers on Mac (e.g.
`http://gcc.gnu.org/wiki/GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_). If
you are unsure about how to do this, speak to your system administrator.

Non-root installs
=================

If you do not have root access to the machine you are using, then replace
``/usr/local`` in the following instructions by e.g. ``$HOME/usr``.
In addition, you should never include ``sudo`` in any of the commands.

Automated Installation
======================

.. note:: You only need to follow this section if you do **not** have HDF5 or
          MPI already installed

The easiest way to install these dependencies correctly is to use the
installation script provided with Hyperion. First, make sure you have
downloaded the latest tar file from `here <https://github.com/hyperion-rt/hyperion/downloads>`_, then expand it with::

    tar xvzf hyperion-x.x.x.tar.gz
    cd hyperion-x.x.x

Then, go to the ``deps/fortran`` directory and run the automated install
script provided::

    cd deps/fortran
    python install.py <prefix>

where ``<prefix>`` is the folder in which you want to install the MPI and HDF5
libraries. To avoid conflicting with existing installed versions (that may not
have Fortran support), it is best to install these in a dedicated directory
such as ``/usr/local/hyperion``::

    python install.py /usr/local/hyperion

and the libraries will be installed in the ``lib``, ``include``, etc. directories inside ``/usr/local/hyperion``. Once the installation is complete, the installer will instruct you to add certain commands to your startup files.

.. note:: if you are installing to a location outside your user directory,
          you will need to run the command with ``sudo``, i.e.::

            sudo python install.py <prefix>

Next, open a new terminal and ensure that the following commands::

    which mpif90
    which h5fc

return a path that is inside the installation path you specified, for example::

    $ which mpif90
    /usr/local/hyperion/bin/mpif90
    $ which h5fc
    /usr/local/hyperion/bin/h5fc

If you get ``command not found`` then you have probably not set up your
``$PATH`` correctly.

The installation script has a number of options (e.g. to set the compilers)
that can be seen with::

    python install.py --help

If the installation fails, a log will be posted to the `Pastebin <http://pastebin.com/>`_ service. Copy the URL and report it either by email or on the Github `Issues <https://www.github.com/hyperion-rt/hyperion/issues>`_.

If the installation succeeds, you can ignore the rest of this document, and move on to the :doc:`python_dependencies`.

Manual Installation: MPI
========================

.. note:: You only need to follow this section if you do **not** have MPI
          already installed.

In order to use the parallel version of the radiation transfer code, you will
need an installation of MPI that supports Fortran. By default, MacOS X ships
with OpenMPI, but the Fortran bindings are not included. In this section, I
have included instructions to install the MPICH2 library with support for
Fortran (though you can in principle use any MPI distribution).

Installation
------------

.. note:: If you encounter any errors at any stage, see the :ref:`mpitrouble` section.

First, download the source for the latest *stable release* of MPICH2 from the
`MPI`_ downloads page. Once downloaded, unpack the file and then go into the
source directory::

    cd mpich2-x.x.x

and configure the installation::

    ./configure --enable-fc --prefix=/usr/local/mpich2

In practice, you will probably want to use a specific fortran compiler, which
you can specify using the ``F77`` and ``FC`` variables as follows::

    ./configure F77=ifort FC=ifort --enable-fc --prefix=/usr/local/mpich2

Once the configure script has successfully run, you can then proceed to build
the MPI library::

    make

If the build is successful, then you can install the library into place using::

    sudo make install

Finally, you will need to add the MPICH2 ``/usr/local/mpich2/bin`` directory to your ``$PATH``.
To check that the installation was successful, type::

    which mpif90

and you should get::

    /usr/local/mpich2/bin/mpif90

If this is not the case, then the installation was unsuccessful.

.. _mpitrouble:

Troubleshooting
---------------

MacOS 10.5 and ifort
^^^^^^^^^^^^^^^^^^^^

If you get the following error when running ./configure::

    configure: error: ****  Incompatible Fortran and C Object File Types!  ****
    F77 Object File Type produced by "ifort  -O2" is : : Mach-O 64-bit object x86_64.
    C  Object File Type produced by "gcc  -O2" is : : Mach-O object i386.

then you are probably using the 64-bit Intel Fortran Compiler on MacOS 10.5.x,
but the 32-bit version of gcc. To fix this, you will need to switch to using
the 32-bit Intel Fortran Compiler. First, clean up the installation so far
with::

    make clean

Then, rerun configure and build using::

    ./configure F77="ifort -m32" FC="ifort -m32" --enable-fc --prefix=/usr/local/mpich2
    make
    sudo make install

Manual Installation: HDF5
=========================

.. note:: You only need to follow this section if you do **not** have HDF5
          already installed.

Installation
------------

.. note:: If you encounter any errors at any stage, see the :ref:`hdftrouble` section.

To compile the Fortran part of the radiation transfer code, you will need the
HDF5 library v1.8.5 or later, with support for Fortran enabled. While package
managers such as Fink and MacPorts include HDF5, they often do not include the
Fortran bindings. Therefore, it is best to install the HDF5 library manually
from source.

To start with, download the source code from the `HDF5 downloads`_ page, then
go into the source code directory::

    cd hdf5-x.x.x

and configure the installation::

    ./configure --enable-fortran --enable-hl --prefix=/usr/local/hdf5_fortran

In practice, you will probably want to use a specific fortran compiler, which
you can specify using the ``FC`` variable as follows::

    ./configure --enable-fortran --enable-hl --prefix=/usr/local/hdf5_fortran FC=ifort

Once the configure script has successfully run, you can then proceed to build
the HDF5 library::

    make

If the build is successful, then you can install the library into place using::

    sudo make install

Finally, you will need to add the HDF5 ``/usr/local/hdf5_fortan/bin`` directory to your ``$PATH``.
To check that the installation was successful, type::

    which h5fc

and you should get::

    /usr/local/hdf5_fortran/bin/h5fc

If this is not the case, then the installation was unsuccessful.

.. note:: The reason we install HDF5 in ``hdf5_fortran`` as opposed to simply
          ``hdf5`` is so as not to conflict with a possible installation of
          the library without the Fortran bindings.


.. _hdftrouble:

Troubleshooting
---------------

MacOS 10.5 and ifort
^^^^^^^^^^^^^^^^^^^^

If you get the following error when running make::

    ...
    H5f90proto.h:1211: warning: previous declaration of 'H5_FC_FUNC_' was here
    H5f90proto.h:1216: error: 'H5_FC_FUNC_' declared as function returning a function
    H5f90proto.h:1216: warning: redundant redeclaration of 'H5_FC_FUNC_'
    H5f90proto.h:1213: warning: previous declaration of 'H5_FC_FUNC_' was here
    H5f90proto.h:1218: error: 'H5_FC_FUNC_' declared as function returning a function
    H5f90proto.h:1218: warning: parameter names (without types) in function declaration
    H5f90proto.h:1218: warning: redundant redeclaration of 'H5_FC_FUNC_'
    H5f90proto.h:1216: warning: previous declaration of 'H5_FC_FUNC_' was here
    make[3]: *** [H5f90kit.lo] Error 1
    make[2]: *** [all] Error 2
    make[1]: *** [all-recursive] Error 1
    make: *** [all-recursive] Error 1

then you are probably using the 64-bit Intel Fortran Compiler on MacOS 10.5.x, but the 32-bit version of gcc.
To fix this, you will need to switch to using the 32-bit Intel Fortran
Compiler. First, clean up the installation so far with::

    make clean

Then, rerun configure and build using::

    ./configure --enable-fortran --enable-hl --prefix=/usr/local/hdf5_fortran FC="ifort -m32"
    make
    sudo make install

If this does not work, try cleaning again, and setup the 32-bit ifort using the scripts provided with ifort. For example, if you are using ifort 11.x, you can do::

    make clean
    source /opt/intel/Compiler/11.0/056/bin/ia32/ifortvars_ia32.sh
    ./configure --enable-fortran --enable-hl --prefix=/usr/local/hdf5_fortran FC=ifort
    make
    sudo make install

NAG f95
^^^^^^^

If you get the following error when running make::

    Error: H5fortran_types.f90, line 39: KIND value (8) does not specify a valid representation method
    Errors in declarations, no further processing for H5FORTRAN_TYPES
    [f95 error termination]
    make[3]: *** [H5fortran_types.lo] Error 1
    make[2]: *** [all] Error 2
    make[1]: *** [all-recursive] Error 1
    make: *** [all-recursive] Error 1

you are using the NAG f95 compiler, which by default does not like statements
like ``real(8) :: a``. To fix this, you will need to specify the
``-kind=byte`` option for the f95 compiler. First, clean up the installation
so far with::

    make clean

Then, rerun configure and build using::

    ./configure --enable-fortran --enable-hl --prefix=/usr/local/hdf5_fortan FC="ifort -kind=byte"
    make
    sudo make install


