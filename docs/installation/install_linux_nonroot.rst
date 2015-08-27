Installing dependencies on Linux without root access
----------------------------------------------------

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

Note that in order to proceed you will need a Fortran compiler such as ``gfortran`` or ``ifort`` (there is no *easy* way to install a Fortran compiler without root access).

To install `HDF5 <http://www.hdfgroup.org/HDF5/>`_ and `MPICH2 <http://www.mpich.org/>`_, which are required by Hyperion, you can do::

    wget https://raw.githubusercontent.com/hyperion-rt/hyperion/master/deps/fortran/install.py
    python install.py $HOME/usr
    
Then, make sure that $HOME/usr/bin is in your ``$PATH``. If you use bash for example, you would do this with::

    export PATH=$HOME/usr/bin:$PATH
  
Next, open a new terminal and ensure that the following commands::

    which mpif90
    which h5fc

return valid paths inside ``$HOME/usr/bin``::

    $ which mpif90
    /home/tom/bin/mpif90
    $ which h5fc
    /home/tom/bin/h5fc

If you get ``command not found`` then you have probably not set up your
``$PATH`` correctly.

The installation script has a number of options (e.g. to set the compilers)
that can be seen with::

    python install.py --help

If the installation fails, a log will be posted to the `Pastebin <http://pastebin.com/>`_ service. Copy the URL and report it either by email or on the Github `Issues <https://www.github.com/hyperion-rt/hyperion/issues>`_.

Python dependencies
^^^^^^^^^^^^^^^^^^^

Download and install the free `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_, which is a scientific linux distribution that contains all the dependencies required by Hyperion. By default, this will be installed to your home directory and will not require root access.

Once you have installed it, you can check that you are indeed using the Python version from Anaconda by typing::

    python --version
    
and you should see something like::

    Python 2.7.10 :: Continuum Analytics, Inc

Hyperion
^^^^^^^^

You are now ready to install Hyperion. Proceed to the :ref:`Hyperion installation instructions <hyperion_install>`!
