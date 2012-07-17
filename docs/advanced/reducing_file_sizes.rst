Reducing file sizes
===================

By default, and to be safe, links between HDF5 files are not used, since they
would potentially break if the files being referred to were moved or deleted.
However, this means that in a few cases, data is duplicated. There are two
places where this occurs, described below.

Linking when calling Model.write
--------------------------------

When writing out a model, the dust properties (and optionally the specific
energy from a previous run) are copied into the input HDF5 file. To use links
instead, you can use::

    m = Model(...)
    ...
    m.write(..., copy=False)

An additional argument, ``absolute_paths``, can be used to specify whether to
use links that are relative (to the input file) or absolute. The default is to
use relative links, but you can use absolute links with::

    m = Model(...)
    ...
    m.write(..., copy=False, absolute_paths=True)

Linking to input files in output files
--------------------------------------

Once the Fortran code has computed an output HDF5 file, it will copy the
contents of the input file into an ``Input/`` group in the output file, as
certain information such as grid geometry can be useful in post-processing. To
link to the input file instead, you can use the following method when setting
up the model in your Python script::

    m = Model(...)
    ...
    m.set_copy_input(False)

In the above example, the Fortran code will now link to the input instead of
copying it. In this case, the path used is the same as the path to the input
file specifying when calling the Hyperion Fortran code.
