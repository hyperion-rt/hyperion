.. _specific:

============================
System-specific instructions
============================

Harvard
=======

Odyssey
-------

You can load up most of the required dependencies for the Fortran and Python code using::

    module load math/hdf5-1.8.4_openmpi-1.3.3_intel-11.1.046
    module load hpc/python-2.6.2

You then will need to install the ATpy module. Since you do not have root access, you will need to install it locally. To do this, create a file called ``.pydistutils.cfg`` in your home directory, and put the following contents inside::

    [install]
    install_lib = ~/usr/python/$py_version_short/site-packages
    install_scripts = ~/usr/bin
    install_data = ~/usr/share
    
Then, add::

    export PYTHONPATH=$HOME/usr/python/2.6/site-packages/
    
to your ``.bashrc`` file.

Then download and install `ATpy <http://atpy.sourceforge.net/>`_ using::

    tar xvzf ATpy-X-X.X.tar.gz
    cd ATpy-X.X.X/
    python setup.py install

You can now proceed to install the Python and Fortran components of the code.

Note that since you do not have access to a display on Odyssey, you will need to make sure matplotlib is in non-interactive mode before importing the ``hyperion`` Python module::

     import matplotlib
     matplotlib.use('Agg')
     import hyperion
    
If you do not do this, you will get the following error::

    RuntimeError: could not create GdkCursor object
