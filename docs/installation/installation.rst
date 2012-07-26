============
Installation
============

.. important:: This section contains information on setting up the dependencies for
          Hyperion as well as Hyperion itself. If you have any issues with the
          installation of any of the dependencies or Hyperion, please first
          talk to your system administrator to see if they can help you get
          set up!

Dependencies
============

First, you will need to install several dependencies. Please follow the
instructions at the following pages to make sure that you have all the
dependencies installed.

.. toctree::
   :maxdepth: 1

   fortran_dependencies.rst
   python_dependencies.rst

.. note:: For instructions for specific computer clusters, see the :ref:`specific` instead, then proceed to the instructions for installing Hyperion below.

.. _hyperion_install:

Hyperion
========

Download the latest tar file from `here <https://github.com/hyperion-rt/hyperion/downloads>`_, then expand it with::

    tar xvzf hyperion-x.x.x.tar.gz
    cd hyperion-x.x.x

Python module
-------------

Install the Python module with::

    python setup.py install

or::

    python setup.py install --user

if you do not have root access. Check that the module installed correctly::

    $ python
    Python 2.7.2 (default, Jan 31 2012, 22:38:06)
    [GCC 4.2.1 (Apple Inc. build 5646)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import hyperion
    >>>

and also try typing::

    $ hyperion

in your shell. If you get ``command not found``, you need to ensure that the
scripts installed by Python are in your ``$PATH``. If you do not know where
these are located, check the last line of the install command above, which
should contain something like this::

    changing mode of /Users/tom/Library/Python/2.7/bin/hyperion to 755

The path listed (excluding ``hyperion`` at the end) should be in your
``$PATH``.

Fortran binaries
----------------

Compile the Fortran code with::

    ./configure
    make
    make install

By default, the binaries will be written to ``/usr/local/bin`` (which will
require you to use ``sudo`` for the last command), but you can change this
using the ``--prefix`` option to configure, for example::

    ./configure --prefix=/usr/local/hyperion

or::

    ./configure --prefix=$HOME/usr

To check that the Fortran binaries are correctly installed, try typing::

    $ hyperion_sph
    Usage: hyperion input_file output_file

If you get::

    $ hyperion_sph
    hyperion_sph: command not found

then something went wrong in the installation, or the directory to which you
installed the binaries is not in your ``$PATH``. Otherwise, you are all set!
