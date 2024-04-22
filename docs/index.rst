.. RT Code documentation documentation master file, created by
   sphinx-quickstart on Fri Feb  5 14:31:14 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Hyperion documentation
======================

Introduction
------------

This is the documentation for `Hyperion <http://www.hyperion-rt.org>`_, a
three-dimensional dust continuum Monte-Carlo radiative transfer code. Models
are set up via Python scripts, and are run using a compiled Fortran code,
optionally making use of the Message Passing Interface (MPI) for parallel
computing.

.. important:: **Before you proceed, please make sure you have read the
               following disclaimers**:

               * The developers cannot guarantee that the code is bug-free,
                 and users should sign up to the `mailing list <https://groups.google.com/forum/?fromgroups#!forum/hyperion-announce>`_ to ensure that
                 they are informed as soon as bugs are identified and fixed.
                 The developers cannot be held responsible for incorrect
                 results, regardless of whether these arise from incorrect
                 usage, a bug in the code, or a mistake in the
                 documentation.

               * Users should read the :doc:`important/important` before using
                 Hyperion. In particular, users are fully responsible for
                 ensuring that parameters such as photon numbers and grid
                 resolution are adequate for the problem being studied.
                 Hyperion will *not* raise errors if these inputs are
                 inadequate.

If your work makes use of Hyperion, please cite:

**Robitaille, 2011**, *HYPERION: an open-source parallelized three-dimensional dust continuum radiative transfer code*, Astronomy & Astrophysics 536 A79 (`ADS <http://adsabs.harvard.edu/abs/2011A%26A...536A..79R>`_, `BibTeX <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2011A%26A...536A..79R&data_type=BIBTEX&db_key=AST&nocookieset=1>`_).

If you need consulting on using Hyperion beyond the current documentation, then
please `contact me <mailto:thomas.robitaille@gmail.com>`_ with details of your
project.

For details of changes between Hyperion versions, see the `changelog <https://github.com/hyperion-rt/hyperion/blob/main/CHANGES.md>`_.

Note on units and constants
---------------------------

All quantities in Hyperion are expressed in the cgs system. Throughout the
documentation, constants are sometimes used in place of values (e.g. ``au``,
``pc``). These can be imported (in Python) using::

    from hyperion.util.constants import *

or, to control which constants are imported::

    from hyperion.util.constants import au, pc, lsun

See :doc:`api/hyperion.util.constants` for more details.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation/installation.rst
   important/important.rst
   setup/setup.rst
   running/running.rst
   postprocessing/postprocessing.rst
   tutorials/index.rst
   dust/dust.rst

Advanced
--------

.. toctree::
   :maxdepth: 1

   advanced/advanced.rst
   api/api.rst
   contributing.rst

Credits
-------

Hyperion is currently being developed by `Thomas Robitaille
<http://www.mpia-hd.mpg.de/~robitaille/>`_.

Interested in contributing fixes or patches to the code or documentation?
Read :doc:`contributing` for more details! If you are interested in
developing new features, `contact me
<mailto:thomas.robitaille@gmail.com>`_ and we can discuss how to coordinate
efforts.

A great thanks to the following users whose help with testing early versions
of Hyperion was invaluable:

* Katharine Johnston
* Nils Lippok
* Stella Offner
* Sarah Ragan
* Andrew Schechtman-Rook
* Amy Stutz
* Barbara Whitney
* Mike Wolff
