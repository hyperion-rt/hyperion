Installing dependencies on MacOS X
----------------------------------

Fortran compiler
^^^^^^^^^^^^^^^^

The first dependency is a Fortran compiler. In addition to commercial
compilers (e.g. ``ifort``, ``pgfortran``, ...), there are a couple of free
ones, the most common of which is ``gfortran``. 
If you don't already have a
compiler installed, you can install ``gfortran`` via package managers on Linux
machines, or from MacPorts or binary installers on Mac (e.g.
`http://gcc.gnu.org/wiki/GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_). If
you are unsure about how to do this, speak to your system administrator.

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

To install `HDF5 <http://www.hdfgroup.org/HDF5/>`_ and `MPICH2 <http://www.mpich.org/>`_, which are required by Hyperion, you can do::

    wget https://raw.githubusercontent.com/hyperion-rt/hyperion/master/deps/fortran/install.py
    sudo python install.py /usr/local

Next, open a new terminal and ensure that the following commands::

    which mpif90
    which h5fc

return the following::

    $ which mpif90
    /usr/local/bin/mpif90
    $ which h5fc
    /usr/local/bin/h5fc

If you get ``command not found`` then you will need to make sure ``/usr/local/bin`` is in your ``$PATH`` (it should be by default though).

The installation script has a number of options (e.g. to set the compilers)
that can be seen with::

    python install.py --help

If the installation fails, a log will be posted to the `Pastebin <http://pastebin.com/>`_ service. Copy the URL and report it either by email or on the Github `Issues <https://www.github.com/hyperion-rt/hyperion/issues>`_.

If you don't have root access, you can do::

    sudo python install.py $HOME/usr

Then be sure to add ``$HOME/usr/bin`` to your ``$PATH``::

    export PATH=$HOME/usr/bin:$PATH

Python dependencies
^^^^^^^^^^^^^^^^^^^

Download and install the free `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_, which is a scientific linux distribution that contains all the dependencies required by Hyperion. By default, this will be installed to your home directory and will not require root access.

Once you have installed it, you can check that you are indeed using the Python
version from Anaconda by typing::

    python --version

and you should see something like::

    Python 2.7.10 :: Continuum Analytics, Inc

Hyperion
^^^^^^^^

You are now ready to install Hyperion. Proceed to the :ref:`Hyperion installation instructions <hyperion_install>`!
