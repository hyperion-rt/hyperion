============
Installation
============

Dependencies
============

First, you will need to install several dependencies. Please follow the
instructions at the following pages to make sure that you have all the
dependencies installed:

.. toctree::
   :maxdepth: 1

   fortran_dependencies.rst
   python_dependencies.rst
   
Hyperion
========

Go into the ``hyperion-x.x.x`` directory, and type::

    python setup.py install
    
to install the Python module, and then::

    ./configure
    make
    sudo make install
    
to install the Fortran code. By default, the binaries will be written to ``/usr/local/bin``, but you can change this using the ``--prefix`` option to configure, for example::

    ./configure --prefix=$HOME/usr
    
To check that the Fortran binaries are correctly installed, try typing::

    $ hyperion_sph
    Usage: bin/rt input_file output_file

If you get::

    $ hyperion_sph
    hyperion_sph: command not found

then something went wrong in the installation, or the directory to which you
installed the binaries is not in your ``$PATH``.


    
