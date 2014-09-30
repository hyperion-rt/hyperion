import numpy as np

from astropy.tests.helper import pytest
from numpy.testing import assert_allclose

from ..model import Model
from .test_helpers import random_id, get_test_dust


class TestSpecificEnergyType(object):

    def setup_method(self, method):

        self.model = Model()

        self.model.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        self.model.set_n_photons(initial=1, imaging=0)
        self.model.set_n_initial_iterations(3)

        self.source = self.model.add_point_source()
        self.source.luminosity = 1.
        self.source.temperature = 1000.

        self.dust = get_test_dust()
        self.density = np.array([[[1.]]])
        self.specific_energy = np.array([[[2.]]])

    def test_specific_energy_initial(self, tmpdir):

        input_file = tmpdir.join(random_id()).strpath
        output_file = tmpdir.join(random_id()).strpath

        self.model.add_density_grid(self.density, self.dust, specific_energy=self.specific_energy)
        self.model.set_minimum_specific_energy(0.5)

        self.model.write(input_file)
        model_out = self.model.run(output_file)

        t = model_out.get_physical_grid('specific_energy')
        assert_allclose(t[0][0, 0, 0], 0.5)

    def test_specific_energy_additional(self, tmpdir):

        input_file = tmpdir.join(random_id()).strpath
        output_file = tmpdir.join(random_id()).strpath

        self.model.add_density_grid(self.density, self.dust, specific_energy=self.specific_energy)
        self.model.set_minimum_specific_energy(0.5)
        self.model.set_specific_energy_type('additional')

        self.model.write(input_file)
        model_out = self.model.run(output_file)

        t = model_out.get_physical_grid('specific_energy')
        assert_allclose(t[0][0, 0, 0], 2.08583984422)

    def test_specific_energy_additional_noiter(self, tmpdir):

        input_file = tmpdir.join(random_id()).strpath
        output_file = tmpdir.join(random_id()).strpath

        self.model.add_density_grid(self.density, self.dust, specific_energy=self.specific_energy)
        self.model.set_minimum_specific_energy(0.5)
        self.model.set_specific_energy_type('additional')

        self.model.set_n_initial_iterations(0)
        self.model.set_n_photons(imaging=0)

        self.model.write(input_file)

        log_file = tmpdir.join(random_id()).strpath

        with pytest.raises(SystemExit) as exc:
            self.model.run(output_file, logfile=log_file)
        assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                    'complete'
        assert "Cannot use specific_energy_type='additional' if the number of" in open(log_file).read()
