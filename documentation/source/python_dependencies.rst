.. _pythondep:

=========================
Python code dependencies
=========================

Overview
========

The packages required for the python code are:

* `NumPy <http://www.scipy.org/>`_
* `SciPy <http://www.scipy.org/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_
* `h5py <http://h5py.alfven.org/>`_

How you install these depends on your operating system, whether you are an existing Python user, and whether you use package managers. To find out whether any of these are already installed, start up a  Python prompt by typing ``python`` on the command line, then try the following commands::

    import numpy
    import scipy
    import matplotlib
    import h5py
    
If you see this::

    >>> import numpy
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: No module named numpy
    >>>

then the module is not installed. If you see this

    >>> import numpy
    >>> 

then the module is already installed.

.. note:: If you do not have root access to the machine you are using, then follow the instructions at `astropython <http://www.astropython.org/tutorial/2010/1/User-rootsudo-free-installation-of-Python-modules>`_ to install the modules into your home directory. In addition, you should not ever include ``sudo`` in any of the commands.

Python newbies
==============

If you want a full Python installation that includes numpy, scipy, matplotlib, h5py, and many other scientific packages, and are a member of an academic institution, but want the easiest solution out there, you can download the `Enthought Python Distribution <http://www.enthought.com/products/edudownload.php>`_ (EPD) for free. EPD includes all the required packages listed above

If you do this, you should be all set, and you can ignore the more detailed installation instructions below.

.. note:: h5py was added to EPD in late 2010, so if you have an older version
          of EPD already installed, you may need to upgrade)

MacOS X
=======

System/python.org Python
------------------------

If you do not want to install EPD, the easiest way to install the three first
packages is to download and install the MacOS X ``dmg`` files for NumPy,
SciPy, and Matplotlib. Use the links at the top of this section to get the
latest dmg files from the different websites. You can of course also install
these from source, but this is beyond the scope of this documentation.

.. note:: If you get an error saying *x can't be installed on this disk. x
          requires Python 2.6 from www.python.org to install*, then this means
          you are probably just using the system Python installation. Go to
          `www.python.org <www.python.org>`_ and download the 2.6.5 or 2.6.6
          version of Python, install, and try installing the packages again.

Check that the packages import correctly::

    $ python
    Python 2.6.1 (r261:67515, Feb 11 2010, 00:51:29) 
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy
    >>> import scipy
    >>> import matplotlib
    >>> 

If any of the packages are incorrectly installed, they will not import cleanly
as above. Once the three main packages are installed, you will need to install
h5py. First, you will need to install the HDF5 library. Note that for the
Fortran code, you also need to install the HDF5 library, but here we need to
create a clean installation without the fortran bindings, or else h5py will
not install properly. Make sure that you perform the following installation in
a different directory from before, to avoid overwriting any files.

To install the plain HDF5 library download the source code from `here
<http://www.hdfgroup.org/ftp/HDF5/current/src/>`_, then go into the source
code directory::

    cd hdf5-1.8.5

and carry out the installation::

    ./configure --prefix=/usr/local/hdf5
    make
    sudo make install

Now, download the latest ``tar.gz`` package from the ``h5py`` `website <http://code.google.com/p/h5py/>`_, and do::

    tar xvzf h5py-x.x.x.tar.gz
    cd h5py-x.x.x
    python setup.py build --api=18 --hdf5=/usr/local/hdf5
    python setup.py install
    
Now, go back to your home directory, and check that ``h5py`` imports cleanly::

    $ python
    Python 2.6.1 (r261:67515, Feb 11 2010, 00:51:29) 
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import h5py
    >>> 
    
    
MacPorts Python
---------------

If you use the Python installation provided by MacPorts, follow the
instructions in this section. To find out if you are using the MacPorts Python
version, type::

    which python
    
If the result is::

    $ which python
    /opt/local/bin/python

you are probably using the MacPorts version. In this case, you can try and
install all of the modules via MacPorts. First, make sure your package list is
up to date::

    sudo port selfupdate

Then do::

    sudo port install py26-numpy py26-scipy py26-matplotlib py26-h5py

If you would prefer to use Python 2.5, replace ``py26`` by ``py25``



