Using grids from HDF5 files
===========================

In the case of large grids, it may in some cases compute once and for all an
HDF5 file containing the grid geometry and quantities, and avoid reading it
when setting up the rest of the model. This can be done using
:meth:`~hyperion.model.Model.use_grid_from_file`, for example::

    from hyperion.model import Model
    m = Model()
    ...
    m.use_grid_from_file('the_grid.hdf5', dust=['kmh.hdf5'])
    
The first argument should be the name of the HDF5 file, and the ``dust``
argument should be a list of dust files to use, one for each dust index in the
grid on disk. Once this has been done, you should not set the grid geometry nor call
:meth:`~hyperion.model.Model.add_density_grid`.

The HDF5 file should be in the format that is written by the ``write`` methods
on grid objects. You can do for example::

    from hyperion.grid import CartesianGrid
    g = CartesianGrid(...)
    # optionally add physical quantities here

then write it out with::

    import h5py
    f = h5py.File('the_grid.hdf5', 'w')
    g.write(f)
    f.close()
