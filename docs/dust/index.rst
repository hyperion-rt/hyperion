Library of dust models for use with Hyperion
============================================

While you are encouraged to prepare your own dust properties file based on
the most appropriate dust model for your problem (see
:doc:`../setup/setup_dust`), we provide a library of dust models that you can
select from. Since dust files tend to be large, the final files are not
directly made available - instead, you should download `this <link>`_ file,
then expand it, and run the ``make_all.py`` script to produce the final dust
files::

    tar xvzf hyperion-dust-0.1.0.tar.gz
    cd hyperion-dust-0.1.0
    make all

The dust files will be generated in the ``dust_files`` directory. If Hyperion
is updated, it may be necessary to re-run this script to update the dust
files to the latest format. If this is necessary, then this will be made
clear in the release announcement for Hyperion.

 .. admonition:: Disclaimer

    The choice of which dust model to use is entirely yours, and you should
    make sure that the dust model that you pick is appropriate for the
    scientific problem you are trying to solve.

The available dust files are:

* ``kmh``: dust properties from `Kim, Martin, and Hendry (1994)`_. The dust
  consists of astronomical silicates, graphite, and carbon and the size
  distribution was determined using a maximum entropy method. This dust type
  is meant to be applicable to the diffuse ISM (for low column densities) in
  the Galaxy.

* More coming soon!

If you develop dust models and would like your dust models to feature as one of the available ones on this page, please let us know!

.. _Kim, Martin, and Hendry (1994): http://dx.doi.org/10.1086/173714
