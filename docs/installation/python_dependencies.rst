.. _pythondep:

=========================
Python code dependencies
=========================

Summary of dependencies
=======================

The packages required for the Python code are:

* `Python <http://www.python.org>`_
* `NumPy <http://www.scipy.org/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_
* `h5py <http://h5py.alfven.org/>`_
* `Astropy <http://www.astropy.org>`_

Overview
========

How you install these depends on your operating system, whether you are an
existing Python user, and whether you use package managers. To find out
whether any of these are already installed, start up a Python prompt by typing
``python`` on the command line, then try the following commands::

    import numpy
    import matplotlib
    import h5py
    import astropy

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

Linux
=====

Numpy, Matplotlib, and h5py should be available in most major Linux package
managers. If Astropy is not available, you can easily install it from source.
To do this, go to the `Astropy homepage <http://www.astropy.org>`_ and
download the latest stable version. Then, expand the tar file and install
using::

    tar xvzf astropy-x.x.tar.gz
    cd astropy-x.x
    python setup.py install

Finally, check that Astropy imports cleanly::

   ~> python
   Python 2.7.2 (default, Aug 19 2011, 20:41:43) [GCC] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
    >>> import astropy
    >>>

MacOS X
=======

MacPorts
--------

If you are installing Python for the first time, we **strongly** recommend the
use of MacPorts to install a full Python distribution. If you would like to do
this, follow `these <http://astrofrog.github.com/macports-python/>`_
instructions to get set up. Once you have your Python distribution installed,
make sure all the dependencies for Hyperion are installed::

    sudo port selfupdate
    sudo port install py27-numpy py27-matplotlib py27-h5py py27-astropy

If this works, you are all set, and you can move on to the actual
:ref:`hyperion_install` installation instructions.

System/python.org Python
------------------------

Numpy and Matplotlib
^^^^^^^^^^^^^^^^^^^^

If you do not want to use MacPorts, the easiest way to install the two first
dependencies is to download and install the MacOS X ``dmg`` files for NumPy
and Matplotlib. Use the links at the top of this section to get the latest
``dmg`` files from the different websites. You can of course also install
these from source, but this is beyond the scope of this documentation.

.. note:: If you get an error saying *x can't be installed on this disk. x
          requires Python 2.7 from www.python.org to install*, then this means
          you are probably just using the system Python installation. Go to
          `www.python.org <http://www.python.org>`_ and download the 2.7.2
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
as above.

h5py
^^^^

Once Numpy and Matplotlib are installed, you will need to install
h5py. First, you will need to install the HDF5 library. Note that for the
Fortran code, you also need to install the HDF5 library, but here we need to
create a clean installation without the fortran bindings, or else h5py will
not install properly. Make sure that you perform the following installation in
a different directory from before, to avoid overwriting any files.

To install the plain HDF5 library download the source code from the latest
`HDF5 downloads <http://www.hdfgroup.org/ftp/HDF5/current/src/>`_ (choose the
hdf5-x.x.x.tar.gz file), then expand the source code::

    tar xvzf hdf5-x.x.x.tar.gz
    cd hdf5-x.x.x

and carry out the installation::

    ./configure --prefix=/usr/local/hdf5
    make
    sudo make install

Now, download the latest ``h5py-x.x.x.tar.gz`` package from the
`h5py website <http://code.google.com/p/h5py/>`_, and do::

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

Astropy
^^^^^^^

Finally, if needed, install Astropy by going to the `Astropy homepage`_ and
downloading the latest stable version. Then, expand the tar file and install
using::

    tar xvzf astropy-x.x.tar.gz
    cd astropy-x.x
    python setup.py install

Finally, check that Astropy imports cleanly::

    $ python
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import astropy
    >>>
