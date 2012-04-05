import pytest
import numpy as np

from .. import Model
from ...grid.amr_grid import AMRGrid, Level, Grid


@pytest.mark.parametrize(('grid_type', 'copy'), [(x, y) for x in ['car', 'sph_pol', 'cyl_pol', 'amr'] for y in [True, False]])
def test_use_quantities_cartesian(grid_type, copy):

    # Set up the initial model
    m = Model()

    # Set up initial geometry
    if grid_type == 'car':
        m.set_cartesian_grid([-1., 1.], [-2., 2.], [-3., 3.])
    elif grid_type == 'cyl':
        m.set_cylindrical_polar_grid([0., 1.], [-1., 1.], [0., 2. * np.pi])
    elif grid_type == 'sph':
        m.set_spherical_polar_grid([0., 1.], [0., np.pi], [0., 2. * np.pi])
    else:
        amr = AMRGrid()
        level = Level()
        grid = Grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 8, 8
        grid.quantities['density'] = np.ones((8, 8, 8))
        level.grids.append(grid)
        amr.levels.append(level)
        m.set_amr_grid(amr)

    # Set up initial density
    if grid_type in ['car', 'cyl', 'sph']:
        m.add_density_grid(np.array([[[1.]]]), 'kmh.hdf5')
    else:
        m.add_density_grid(amr['density'], 'kmh.hdf5')

    # Set up source
    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 1000.

    # Set number of photons
    m.set_n_photons(initial=1000, imaging=1000)

    # Write out and run the model
    m.write('test_%s.rtin' % grid_type, overwrite=True, copy=copy)
    m.run('test_%s.rtout' % grid_type, overwrite=True)

    # Set up a second model that uses the properties from the first
    m2 = Model()

    # Use the geometry from the initial model
    m2.use_geometry('test_%s.rtout' % grid_type)

    # Use the density and specific_energy from the initial model
    m2.use_quantities('test_%s.rtout' % grid_type)

    # Set up source again
    s = m2.add_point_source()
    s.luminosity = 1.
    s.temperature = 1000.

    # Set number of photons
    m2.set_n_photons(initial=1000, imaging=1000)

    # Write out and run to test that the file is coherent
    m2.write('test_%s_2.rtin' % grid_type, overwrite=True, copy=copy)
    m2.run('test_%s_2.rtout' % grid_type, overwrite=True)
