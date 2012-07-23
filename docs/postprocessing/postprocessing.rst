.. _post-processing:

======================
Post-processing models
======================

You've successfully set up and run your model, but now you want to extract the results and do something useful with them. The following sections describe how to retrieve observables (SEDs and images) and quantities (such as density and temperature).

.. toctree::
   :maxdepth: 1

   extracting_observables
   extracting_quantities

The :doc:`../tutorials/index` section contains a number of examples illustrating how to extract and visualize observables and quantities.

.. note:: A convenience script is provided to quickly extract image cubes
          and physical grids for inspection as FITS files::

              # Retrieve all images
              hyperion2fits --images model.rtout

              # Retrieve all physical arrays (only works for regular grids)
              hyperion2fits --physics model.rtout

          Ds9 version 7 or later is required to be able to navigate FITS
          cubes with more than 3 dimensions. This is only meant as a quick
          look tool. To extract properly scaled and sliced information
          from the output file, see the sections above.

