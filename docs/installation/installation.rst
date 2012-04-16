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

Hyperion
========

Go into the ``hyperion-x.x.x`` directory, and type::

    python setup.py install

to install the Python module, and then::

    ./configure
    make
    make install

to install the Fortran code. By default, the binaries will be written to ``/usr/local/bin``, but you can change this using the ``--prefix`` option to configure, for example::

    ./configure --prefix=$HOME/usr

To check that the Fortran binaries are correctly installed, try typing::

    $ hyperion_sph
    Usage: hyperion input_file output_file

If you get::

    $ hyperion_sph
    hyperion_sph: command not found

then something went wrong in the installation, or the directory to which you
installed the binaries is not in your ``$PATH``.
