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
instance. These are described in the following sections (the sections on
images and radiative transfer configuration are the same as for the :class:`~hyperion.model.Model`
class)

.. toctree::
   :maxdepth: 1

   setup_dust
   setup_yso
   setup_images
   setup_conf

Once the model is set up, the user can write it out to the disk for use
with the Fortran radiation transfer code::

    m.write('example.rtin')
    