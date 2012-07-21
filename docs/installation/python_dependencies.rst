.. _pythondep:

=========================
Python code dependencies
=========================

Overview
========

The packages required for the Python code are:

* `NumPy <http://www.scipy.org/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_
* `h5py <http://h5py.alfven.org/>`_
* `ATpy <http://atpy.github.com>`_

How you install these depends on your operating system, whether you are an existing Python user, and whether you use package managers. To find out whether any of these are already installed, start up a  Python prompt by typing ``python`` on the command line, then try the following commands::

    import numpy
    import matplotlib
    import h5py
    import atpy

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

.. note:: If you are installing Hyperion from the git repository, you will
          also need the `Cython <http://www.cython.org>`_ module installed.

MacOS X
=======

MacPorts
--------

If you are installing Python for the first time, we **strongly** recommend the use of MacPorts to install a full Python distribution. If you would like to do this, follow `these <http://astrofrog.github.com/macports-python/>`_ instructions to get set up. Once you have your Python distribution installed, make sure all the dependencies for Hyperion are installed::

    sudo port selfupdate
    sudo port install py27-numpy py27-matplotlib py27-h5py py27-atpy

System/python.org Python
------------------------

If you do not want to use MacPorts, the easiest way to install the three first
dependencies is to download and install the MacOS X ``dmg`` files for NumPy
and Matplotlib. Use the links at the top of this section to get the latest
``dmg`` files from the different websites. You can of course also install
these from source, but this is beyond the scope of this documentation.

.. note:: If you get an error saying *x can't be installed on this disk. x
          requires Python 2.7 from www.python.org to install*, then this means
          you are probably just using the system Python installation. Go to
          `www.python.org <www.python.org>`_ and download the 2.7.2
          version of Python, install, and try installing the packages again.

Check that the packages import correctly::

    $ python
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy
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

    cd hdf5-x.x.x

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
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import h5py
    >>>

Finally, if needed, install ATpy by going `here <http://atpy.github.com/>`_ and clicking on *Download latest stable version*. Then, expand the tar file and install using::

    tar xvzf ATpy-x.x.x.tar.gz
    cd ATpy-x.x.x
    python setup.py install

Finally, check that ATpy imports cleanly::

    $ python
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import atpy
    >>>
