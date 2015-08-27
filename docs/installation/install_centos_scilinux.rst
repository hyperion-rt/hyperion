Installing dependencies on CentOS and Scientific Linux (requires root)
----------------------------------------------------------------------

.. note:: The HDF5 and Astropy packages are not available for these
          distributions, so a few additional steps are needed.

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

The Fortran Hyperion code requires a Fortran compiler, `HDF5 <http://www.hdfgroup.org/HDF5/>`_ and `MPICH2 <http://www.mpich.org/>`_.

You can install some of these dependencies with::

    sudo yum install mpich-devel gcc-gfortran libgfortran-static
    export PATH=/usr/lib64/mpich/bin:$PATH

Then you can install HDF5 with::

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

If instead you use the `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_ you can install the Python dependencies with::

    conda install numpy matplotlib h5py astropy
    
Hyperion
^^^^^^^^

You are now ready to install Hyperion. Proceed to the :ref:`Hyperion installation instructions <hyperion_install>`!