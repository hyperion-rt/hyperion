Installing dependencies with Linux package managers
===================================================

Debian/Ubuntu
-------------

You can install the dependencies for the Fortran code with::

    sudo apt-get install libmpich2-dev libhdf5-serial-dev

and the dependencies for the Python code with::

    sudo apt-get install python-numpy python-dev python-astropy python-h5py python-matplotlib

Once you have installed these, you can proceed to the :ref:`Hyperion installation instructions <hyperion_install>`

Red Hat/Fedora/Scientific Linux/CentOS
--------------------------------------

You can install the dependencies for the Fortran code with::

    sudo yum install hdf5-static mpich-devel libgfortran-static
    export PATH=/usr/lib64/mpich/bin:$PATH

and the dependencies for the Python code with::

    sudo yum install gcc-c++ numpy h5py python-matplotlib python-astropy
