Installing dependencies with Linux package managers
===================================================

Debian/Ubuntu
-------------

You can install the dependencies for the Fortran code with::

    sudo apt-get install libmpich2-dev libhdf5-serial-dev

and the dependencies for the Python code with::

    sudo apt-get install python-numpy python-dev python-astropy python-h5py python-matplotlib

Once you have installed these, you can proceed to the :ref:`Hyperion installation instructions <hyperion_install>`

Fedora
------

You can install the dependencies for the Fortran code with::

    sudo yum install hdf5-static mpich-devel libgfortran-static
    export PATH=/usr/lib64/mpich/bin:$PATH

and the dependencies for the Python code with::

    sudo yum install gcc-c++ numpy h5py python-matplotlib python-astropy

CentOS and Scientific Linux
---------------------------

.. note:: The HDF5 and Astropy packages are not available for these
          distributions, so a few additional steps are needed.

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

You can install some of the dependencies for the Fortran code with::

    sudo yum install mpich-devel gcc-gfortran libgfortran-static
    export PATH=/usr/lib64/mpich/bin:$PATH

Now we can install HDF5::

    wget https://raw.githubusercontent.com/hyperion-rt/hyperion/master/deps/fortran/install.py
    sudo python install.py /usr/local --only-hdf5

You can ignore the message about setting the ``PATH``, since ``/usr/local/bin`` should already be in your ``PATH`` by default.

Python dependencies
^^^^^^^^^^^^^^^^^^^

If you use the system Python installation, you can install some of the
dependencies for the Python code with::

    sudo yum install gcc-c++ python-devel numpy python-matplotlib Cython

then you can install h5py and Astropy::

    sudo easy_install pip
    pip install six --user --upgrade
    pip install h5py astropy --user

If instead you use the [Anaconda Python Distribution](https://store.continuum.io/cshop/anaconda/) you can install the Python dependencies with::

    conda install numpy matplotlib h5py astropy

