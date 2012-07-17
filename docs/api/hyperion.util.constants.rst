=================================
hyperion.util.constants
=================================

The ``hyperion.util.constants`` module contains a number of useful physical
constants that can be used when setting up or postprocessing models. All
constants are defined in the cgs system. Since all attributes for objects
(sources, images, etc.) required parameters to be specified in the cgs system,
one can therefore do::

    from hyperion.util.constants import au, pc
    e = m.add_power_law_envelope()
    e.rmin = 0.1 * au
    e.rmax = pc

which is equivalent to writing::

    e = m.add_power_law_envelope()
    e.rmin = 1.49598e12
    e.rmax = 3.08568025e18

but the former is more readable. The available constants are described below:

.. currentmodule:: hyperion.util.constants

Fundamental constants
---------------------

.. autodata:: h
.. autodata:: k
.. autodata:: c
.. autodata:: G
.. autodata:: sigma
.. autodata:: m_h

Solar constants
---------------

.. autodata:: lsun
.. autodata:: msun
.. autodata:: rsun
.. autodata:: tsun

Common Astronomical constants
-----------------------------

.. autodata:: au
.. autodata:: pc
.. autodata:: kpc
.. autodata:: year
