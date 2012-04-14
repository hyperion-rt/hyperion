Introduction to the Python interface
====================================

This section presents an introduction to the Python interface to the
radiation transfer code. The Python module used to interface with the code
is named ``hyperion``, and can be imported using::

    import hyperion

However, most of the useful functionality is stored in sub-modules. For
example, the ``hyperion.model`` module contains classes to easily set up
models, while the ``hyperion.util.constants`` contains useful constants in
cgs units::

    from hyperion.util.constants import pi, lsun, rsun, au, pc