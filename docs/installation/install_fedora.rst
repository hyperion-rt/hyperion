Installing dependencies on Fedora (requires root)
-------------------------------------------------

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

The Fortran Hyperion code requires a Fortran compiler, `HDF5 <http://www.hdfgroup.org/HDF5/>`_ and `MPICH2 <http://www.mpich.org/>`_.

You can install these dependencies with::

    sudo yum install hdf5-static mpich-devel libgfortran-static
    export PATH=/usr/lib64/mpich/bin:$PATH

Python dependencies
^^^^^^^^^^^^^^^^^^^

If you use the system Python installation, you can install the dependencies for
the Python code with::

    sudo yum install gcc-c++ numpy h5py python-matplotlib python-astropy

If instead you use the `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_ you can install the Python dependencies with::

    conda install numpy matplotlib h5py astropy

Hyperion
^^^^^^^^

You are now ready to install Hyperion. Proceed to the :ref:`Hyperion installation instructions <hyperion_install>`!