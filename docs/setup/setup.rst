.. _setup:

Setting up models
-----------------

The easiest way to set up models is via the Hyperion Python package.

You should start by choosing the type of model you want to set up. At the
moment, you can either set up an arbitrary model (which allows you to use an
arbitrary grid and density structure), or an analytical YSO model, which is
specifically models with fixed density structures such as disks, envelopes,
bipolar cavities, and defined on a spherical or cylindrical polar grid. Other
kinds of convenience models may be added in future.

Once you have decided on the type of model, you will need to set up the grid,
sources, dust properties, density structure, image and SED parameters, and
choose the settings for the radiative transfer algorithm.

The following pages give instructions on setting up the two main kinds of
models:

.. toctree::
   :maxdepth: 2

   model.rst
   analytical_yso_model.rst

.. note:: it is possible to write model input files yourself directly in HDF5
          and bypass the Python library entirely (but this is reserved for
          advanced users!). See :doc:`../advanced/advanced` for more
          information.