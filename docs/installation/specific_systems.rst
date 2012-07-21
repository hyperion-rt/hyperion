:orphan:

.. _specific:

============================
System-specific instructions
============================

Harvard/Odyssey
===============

Installing
----------

You can load up most of the required dependencies for the Fortran and Python code using::

    module load math/hdf5-1.8.4_openmpi-1.3.3_intel-11.1.046
    module load hpc/python-2.6.2

Then download and install ATpy by going `here <http://atpy.github.com/>`_ and
clicking on *Download latest stable version*. Then, expand the tar file and
install using::

    tar xvzf ATpy-x.x.x.tar.gz
    cd ATpy-x.x.x
    python setup.py install

You can now proceed to install the Python and Fortran components of the code.

Running
-------

Note that since you do not have access to a display on Odyssey, you will need to make sure matplotlib is in non-interactive mode before importing the ``hyperion`` Python module::

     import matplotlib
     matplotlib.use('Agg')
     import hyperion

If you do not do this, you will get the following error::

    RuntimeError: could not create GdkCursor object
