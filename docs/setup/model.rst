.. _model:

================
Arbitrary Models
================

.. note:: The current document only shows example use of some methods, and
          does not discuss all the options available. To see these, do not
          hesitate to use the ``help`` command, for example ``help
          m.write`` will return more detailed instructions on using the
          ``write`` method.

To create a general model, you will need to first import the Model class
from the Python Hyperion module::

    from hyperion.model import Model

it is then easy to set up a generic model using::

    m = Model()

The model can then be set up using methods of the :class:`~hyperion.model.Model` instance. These
are described in the following sections.

.. toctree::
   :maxdepth: 1

   setup_dust
   setup_grid
   setup_sources
   setup_images
   setup_conf

Once the model is set up, you can write it out to the disk for use
with the Fortran radiation transfer code::

    m.write('example.rtin')

See :meth:`~hyperion.model.Model.write` for information about the
available options.
