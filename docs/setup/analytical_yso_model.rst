.. _analyticalysomodel:

=====================
Analytical YSO Models
=====================

The :ref:`model` class should allow users to set up arbitrary problems.
However, the Python module also provides classes that build on top of
:class:`~hyperion.model.Model` that make it easy to specify certain kinds of problems. These
classes support all the methods available for the :class:`~hyperion.model.Model` class, and define
new ones. The :class:`~hyperion.model.AnalyticalYSOModel` makes it easy to set up sources with
flared disks, and rotationally flattened envelopes, optionally with bipolar
cavities. To use this class, you will first need to import it::

    from hyperion.model import AnalyticalYSOModel

it is then easy to set up such a model using::

    m = AnalyticalYSOModel()

The model can then be set up using methods of the :class:`~hyperion.model.AnalyticalYSOModel`
instance. This is described in more detail in the following section:

.. toctree::
   :maxdepth: 1

   setup_yso

To set up the dust, images, and configuration, see the :doc:`setup_dust`,
:doc:`setup_images`, and :doc:`setup_conf` sections of the :doc:`model`
description.

Once the model is set up, you can write it out to the disk for use
with the Fortran radiation transfer code::

    m.write('example.rtin')

See :meth:`~hyperion.model.AnalyticalYSOModel.write` for information about the
available options.

.. note:: One of the available options is ``merge_if_possible=``, which if set
          to ``True`` will merge the various density components into a single
          one if the dust types match. This allows the code to run faster, but
          on the other hand means that if tracking e.g. photon origin, the
          separate origin of components that have been merged will be lost.
          This option is enabled by default, but you may want to disable it.
    