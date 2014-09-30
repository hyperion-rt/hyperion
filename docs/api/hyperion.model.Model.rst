=================================
hyperion.model.Model
=================================

.. currentmodule:: hyperion.model

.. autoclass:: Model

   .. rubric:: Adding sources

   .. autosummary::

      ~Model.add_source
      ~Model.add_point_source
      ~Model.add_spherical_source
      ~Model.add_external_spherical_source
      ~Model.add_external_box_source
      ~Model.add_map_source
      ~Model.add_plane_parallel_source

   .. rubric:: Setting the grid

   .. autosummary::

      ~Model.set_grid
      ~Model.set_cartesian_grid
      ~Model.set_cylindrical_polar_grid
      ~Model.set_spherical_polar_grid
      ~Model.set_octree_grid
      ~Model.set_amr_grid

   .. rubric:: Setting quantities

   .. autosummary::

      ~Model.add_density_grid

   .. rubric:: Images/SEDs

   .. autosummary::

      ~Model.add_peeled_images
      ~Model.add_binned_images

   .. rubric:: Configuration

   .. autosummary::

      ~Model.set_seed
      ~Model.set_propagation_check_frequency
      ~Model.set_monochromatic
      ~Model.set_minimum_temperature
      ~Model.set_minimum_specific_energy
      ~Model.set_specific_energy_type
      ~Model.set_n_initial_iterations
      ~Model.set_raytracing
      ~Model.set_max_interactions
      ~Model.set_max_reabsorptions
      ~Model.set_pda
      ~Model.set_mrw
      ~Model.set_convergence
      ~Model.set_kill_on_absorb
      ~Model.set_forced_first_scattering
      ~Model.set_output_bytes
      ~Model.set_sample_sources_evenly
      ~Model.set_enforce_energy_range
      ~Model.set_copy_input

   .. rubric:: Running

   .. autosummary::

      ~Model.write
      ~Model.run

   .. rubric:: Using grids from files

   .. autosummary::

      ~Model.use_grid_from_file

   .. rubric:: Re-using previous models

   .. autosummary::

      ~Model.read
      ~Model.use_geometry
      ~Model.use_quantities
      ~Model.use_sources
      ~Model.use_image_config
      ~Model.use_run_config
      ~Model.use_output_config

   .. rubric:: Methods (detail)

   .. automethod:: Model.add_source
   .. automethod:: Model.add_point_source
   .. automethod:: Model.add_spherical_source
   .. automethod:: Model.add_external_spherical_source
   .. automethod:: Model.add_external_box_source
   .. automethod:: Model.add_map_source
   .. automethod:: Model.add_plane_parallel_source

   .. automethod:: Model.set_grid
   .. automethod:: Model.set_cartesian_grid
   .. automethod:: Model.set_cylindrical_polar_grid
   .. automethod:: Model.set_spherical_polar_grid
   .. automethod:: Model.set_octree_grid
   .. automethod:: Model.set_amr_grid

   .. automethod:: Model.add_density_grid

   .. automethod:: Model.add_peeled_images
   .. automethod:: Model.add_binned_images

   .. automethod:: Model.set_seed
   .. automethod:: Model.set_propagation_check_frequency
   .. automethod:: Model.set_monochromatic
   .. automethod:: Model.set_minimum_temperature
   .. automethod:: Model.set_minimum_specific_energy
   .. automethod:: Model.set_specific_energy_type
   .. automethod:: Model.set_n_initial_iterations
   .. automethod:: Model.set_raytracing
   .. automethod:: Model.set_max_interactions
   .. automethod:: Model.set_max_reabsorptions
   .. automethod:: Model.set_pda
   .. automethod:: Model.set_mrw
   .. automethod:: Model.set_convergence
   .. automethod:: Model.set_kill_on_absorb
   .. automethod:: Model.set_forced_first_scattering
   .. automethod:: Model.set_output_bytes
   .. automethod:: Model.set_sample_sources_evenly
   .. automethod:: Model.set_enforce_energy_range
   .. automethod:: Model.set_copy_input

   .. automethod:: Model.write
   .. automethod:: Model.run

   .. automethod:: Model.use_grid_from_file

   .. automethod:: Model.read
   .. automethod:: Model.use_geometry
   .. automethod:: Model.use_quantities
   .. automethod:: Model.use_sources
   .. automethod:: Model.use_image_config
   .. automethod:: Model.use_run_config
   .. automethod:: Model.use_output_config
