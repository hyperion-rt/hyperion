===========================
High-Level Python Interface
===========================

The main core of the radiation transfer code is written in fortran. This part of the code requires a single HDF5 file which describes the geometry and physics of the problem. Creating this file from scratch can be tedious, so a high-level interface is provided, written in Python.

To import this high level interface, simply use::

    import hyperion

