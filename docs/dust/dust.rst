Library of dust models
======================

While you are encouraged to prepare your own dust properties file based on the
most appropriate dust model for your problem (see :doc:`../setup/setup_dust`),
we provide a library of dust models that you can select from. Since dust files
tend to be large, the final files are not directly made available - instead,
you should download `this`_ file, then expand it, and run the ``setup.py``
script as shown to produce the final dust files::

    tar xvzf hyperion-dust-0.1.0.tar.gz
    cd hyperion-dust-0.1.0
    python setup.py build_dust

The dust files will be generated in the ``dust_files`` directory. If Hyperion
is updated, it may be necessary to re-run this script to update the dust
files to the latest format. If this is necessary, then this will be made
clear in the release announcement for Hyperion.

 .. admonition:: Disclaimer

    The choice of which dust model to use is entirely yours, and you should
    make sure that the dust model that you pick is appropriate for the
    scientific problem you are trying to solve. In addition, as for the
    radiative transfer code itself, we cannot guarantee that the dust files are
    bug-free - you should sign up to the `mailing list
    <https://groups.google.com/forum/?fromgroups#!forum/hyperion-announce>`_ to
    ensure that you are informed as soon as bugs are identified and fixed.


The available dust files at this time are:

.. toctree::
   :maxdepth: 1

   kmh_hg
   kmh
   d03

If you develop dust models and would like your dust models to feature as one of
the available ones on this page, please let us know!

For advanced users, we provide code and documentation to compute your own dust
models:

.. toctree::
   :maxdepth: 1

   bhmie

.. _this: http://pypi.python.org/packages/source/h/hyperion-dust/hyperion-dust-0.1.0.tar.gz
