.. _MPI: http://www.mpich.org/downloads/
.. _HDF5 downloads: http://www.hdfgroup.org/ftp/HDF5/current/src/


Installing dependencies the hard way
------------------------------------

Summary
^^^^^^^

The packages required for the Fortran code are:

* A recent Fortran compiler. The following compilers are known to work:

  * gfortran 4.3 and later
  * ifort 11 and later
  * pgfortran 11 and above

* `HDF5 <http://www.hdfgroup.org/HDF5/>`_ 1.8.5 or later with the Fortran bindings

* An MPI installation (e.g. `MPICH2 <http://www.mpich.org/>`_ or `OpenMPI <http://www.open-mpi.org/>`_) with the Fortran bindings

Note that often, default installations of HDF5 and MPI packages do not include support for Fortran - this has to be explicitly enabled as described below.

The packages required for the Python code are:

* `Python <http://www.python.org>`_
* `NumPy <http://www.scipy.org/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_
* `h5py <http://h5py.alfven.org/>`_
* `Astropy <http://www.astropy.org>`_

MPI
^^^

.. note:: You only need to follow this section if you do **not** have MPI
          already installed.

In order to use the parallel version of the radiation transfer code, you will
need an installation of MPI that supports Fortran. By default, MacOS X ships
with OpenMPI, but the Fortran bindings are not included. In this section, I
have included instructions to install the MPICH2 library with support for
Fortran (though you can in principle use any MPI distribution).

.. note:: If you encounter any errors at any stage, see the :ref:`mpitrouble` section.

First, download the source for the latest *stable release* of MPICH2 from the
`MPI`_ downloads page. Once downloaded, unpack the file and then go into the
source directory::

    cd mpich2-x.x.x

and configure the installation::

    ./configure --enable-fc --prefix=/usr/local

In practice, you will probably want to use a specific fortran compiler, which
you can specify using the ``F77`` and ``FC`` variables as follows::

    ./configure F77=ifort FC=ifort --enable-fc --prefix=/usr/local

Once the configure script has successfully run, you can then proceed to build
the MPI library::

    make

If the build is successful, then you can install the library into place using::

    sudo make install

Finally, if this is not alrady the case, you will need to add the MPICH2
``/usr/local/bin`` directory to your ``$PATH``. To check that the installation
was successful, type::

    which mpif90

and you should get::

    /usr/local/bin/mpif90

If this is not the case, then the installation was unsuccessful.

.. _mpitrouble:

Troubleshooting
"""""""""""""""

MacOS 10.5
**********

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

    ./configure F77="ifort -m32" FC="ifort -m32" --enable-fc --prefix=/usr/local
    make
    sudo make install

HDF5
^^^^

.. note:: You only need to follow this section if you do **not** have HDF5
          already installed.

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
"""""""""""""""

MacOS 10.5 and ifort
********************

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
*******

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

Python
^^^^^^

How you install the Python dependencies depends on your operating system,
whether you are an existing Python user, and whether you use package managers.
To find out whether any of these are already installed, start up a Python
prompt by typing ``python`` on the command line, then try the following
commands::

    import numpy
    import matplotlib
    import h5py
    import astropy

If you see this::

    >>> import numpy
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: No module named numpy
    >>>

then the module is not installed. If you see this

    >>> import numpy
    >>>

then the module is already installed.

Anaconda
""""""""

If you are not sure how to set up a Scientific Python environment, simply download and install the `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_, which contains all the required Python dependencies for Hyperion.

Linux
"""""

Numpy, Matplotlib, and h5py should be available in most major Linux package
managers. If Astropy is not available, you can easily install it from source with::

    pip install astropy --user

MacOS X
"""""""

If you use MacPorts, you can install all the dependencies for Hyperion with::

    sudo port selfupdate
    sudo port install py27-numpy py27-matplotlib py27-h5py py27-astropy


Hyperion
^^^^^^^^

You are now ready to install Hyperion. Proceed to the :ref:`Hyperion installation instructions <hyperion_install>`!