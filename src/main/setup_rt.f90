module setup

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use grid_physics, only : setup_grid_physics
  use grid_geometry, only : setup_grid_geometry
  use sources
  use dust_main
  use type_dust
  use lib_conf
  use binned_images
  use peeled_images
  use settings
  use forced_interaction, only : WR99, BAES16, baes16_xi

  implicit none
  save

  private
  public :: setup_initial
  public :: setup_final_iteration

contains

  subroutine setup_initial(input_handle)

    implicit none

    integer(hid_t),intent(in) :: input_handle
    integer(hid_t) :: g_dust, g_geometry, g_physics, g_sources, g_output
    integer :: physics_io_bytes
    type(version) :: python_version
    integer :: idust
    character(len=10) :: forced_first_interaction_algorithm_str

    if(mp_exists_keyword(input_handle, '/', 'python_version')) then
       call mp_read_keyword(input_handle, '/', 'python_version', python_version%string)
       if(python_version < version('0.8.7')) then
          if(main_process()) call error("setup_initial", "cannot read files made with the Python module before version 0.8.7")
       end if
    else
       if(main_process()) call error("setup_initial", "cannot read files made with the Python module before version 0.8.7")
    end if

    call mp_join()

    call mp_read_keyword(input_handle, '/', 'monochromatic', use_exact_nu)
    if(use_exact_nu) then
        if(mp_exists_keyword(input_handle, '/', 'monochromatic_energy_threshold')) then
            call mp_read_keyword(input_handle, '/', 'monochromatic_energy_threshold', monochromatic_energy_threshold)
        else
            monochromatic_energy_threshold = 1e-10_dp
        end if
    end if
    call mp_read_keyword(input_handle, '/', 'raytracing', use_raytracing)

    call mp_read_keyword(input_handle, '/', 'n_stats', n_stats)

    call mp_read_keyword(input_handle, '/', 'n_inter_max', n_inter_max)
    if(mp_exists_keyword(input_handle, '/', 'n_inter_max_warn')) then
       call mp_read_keyword(input_handle, '/', 'n_inter_max_warn', n_inter_max_warn)
    else
       n_inter_max_warn = .true.
    end if

    call mp_read_keyword(input_handle, '/', 'n_reabs_max', n_reabs_max)
    if(mp_exists_keyword(input_handle, '/', 'n_reabs_max_warn')) then
       call mp_read_keyword(input_handle, '/', 'n_reabs_max_warn', n_reabs_max_warn)
    else
       n_reabs_max_warn = .true.
    end if

    call mp_read_keyword(input_handle, '/', 'pda', use_pda)
    call mp_read_keyword(input_handle, '/', 'mrw', use_mrw)

    if(use_mrw) then
       call mp_read_keyword(input_handle, '/', 'mrw_gamma', mrw_gamma)
       call mp_read_keyword(input_handle, '/', 'n_inter_mrw_max', n_mrw_max)
       if(mp_exists_keyword(input_handle, '/', 'n_inter_mrw_max_warn')) then
          call mp_read_keyword(input_handle, '/', 'n_inter_mrw_max_warn', n_mrw_max_warn)
       else
          n_mrw_max_warn = .true.
       end if
    end if

    call mp_read_keyword(input_handle, '/', 'kill_on_absorb', kill_on_absorb)
    if(mp_exists_keyword(input_handle, '/', 'kill_on_scatter')) then
       call mp_read_keyword(input_handle, '/', 'kill_on_scatter', kill_on_scatter)
    else
       kill_on_scatter = .false.
    end if

    if(mp_exists_keyword(input_handle, '/', 'forced_first_scattering')) then
      call mp_read_keyword(input_handle, '/', 'forced_first_scattering', forced_first_interaction)
      forced_first_interaction_algorithm = WR99
    else
      call mp_read_keyword(input_handle, '/', 'forced_first_interaction', forced_first_interaction)
      call mp_read_keyword(input_handle, '/', 'forced_first_interaction_algorithm', forced_first_interaction_algorithm_str)
      select case(trim(forced_first_interaction_algorithm_str))
      case('wr99')
         forced_first_interaction_algorithm = WR99
      case('baes16')
         forced_first_interaction_algorithm = BAES16
         call mp_read_keyword(input_handle, '/', 'forced_first_interaction_baes16_xi', baes16_xi)
      case default
         call error('setup_initial', 'Unknown forced first interaction algorithm: '//trim(forced_first_interaction_algorithm_str))
      end select
    end if

    if(mp_exists_keyword(input_handle, '/', 'propagation_check_frequency')) then
       call mp_read_keyword(input_handle, '/', 'propagation_check_frequency', propagation_check_frequency)
    else
       propagation_check_frequency = 1.e-3
    end if

    if(mp_exists_keyword(input_handle, '/', 'sample_sources_evenly')) then
       call mp_read_keyword(input_handle, '/', 'sample_sources_evenly', sample_sources_evenly)
    else
       sample_sources_evenly = .false.
    end if

    ! ENERGY RANGE (needs to be set before running the physics set-up)
    if(mp_exists_keyword(input_handle, '/', 'enforce_energy_range')) then
       call mp_read_keyword(input_handle, '/', 'enforce_energy_range', enforce_energy_range)
    else
       enforce_energy_range = .true.
    end if

    ! DUST

    g_dust = mp_open_group(input_handle, '/Dust')
    call setup_dust(g_dust)
    call mp_close_group(g_dust)

    if(n_dust==0) then
       call warn("main", "no dust present, so skipping initial iterations")
       n_initial_iter=0
       if(use_exact_nu) n_last_photons_dust = 0
       if(use_raytracing) n_raytracing_photons_dust = 0
    else
       call mp_read_keyword(input_handle, '/', 'n_initial_iter', n_initial_iter)
       if(n_initial_iter > 0) then
          call mp_read_keyword(input_handle, '/', 'n_initial_photons', n_initial_photons)
          if(n_initial_photons==0) call error("setup_initial", "Number of initial iterations is non-zero, but number of specific_energy photons is zero")
       else
          n_initial_photons = 0
       end if
       if(use_exact_nu) call mp_read_keyword(input_handle, '/', 'n_last_photons_dust', n_last_photons_dust)
       if(use_raytracing) call mp_read_keyword(input_handle, '/', 'n_ray_photons_dust', n_raytracing_photons_dust)
    end if

    ! GRID

    ! Read in specific energy type
    if(mp_exists_keyword(input_handle, '/', 'specific_energy_type')) then
       call mp_read_keyword(input_handle, '/', 'specific_energy_type', specific_energy_type)
    else
       specific_energy_type = 'initial'
    end if

    if (trim(specific_energy_type) == 'additional') then
       if (n_initial_iter == 0) then
          call error("setup_initial", "Cannot use specific_energy_type='additional' if the number of specific energy iterations is 0")
       end if
    else if (trim(specific_energy_type)/='initial') then
       call error("setup_initial", "specific_energy_type should be 'additional' or 'initial'")
    end if

    g_geometry = mp_open_group(input_handle, '/Grid/Geometry')
    call setup_grid_geometry(g_geometry)
    call mp_close_group(g_geometry)

    g_physics = mp_open_group(input_handle, '/Grid/Quantities')
    call setup_grid_physics(g_physics, use_mrw, use_pda)
    call mp_close_group(g_physics)

    call mp_read_keyword(input_handle, '/', 'physics_io_bytes', physics_io_bytes)

    select case(physics_io_bytes)
    case(4)
       physics_io_type = sp
    case(8)
       physics_io_type = dp
    case default
       call error("setup_initial", "unexpected value of physics_io_bytes (should be 4 or 8)")
    end select

    ! FREQUENCIES

    if(use_exact_nu) then
       call mp_table_read_column_auto(input_handle, 'frequencies', 'nu', frequencies)
    end if

    ! SOURCES

    g_sources = mp_open_group(input_handle, '/Sources')
    call setup_sources(g_sources)
    call mp_close_group(g_sources)

    ! If no sources have been set up, give an error if we are not in raytracing only mode
    if(n_sources == 0) then
       if(n_initial_iter > 0) call error("setup_initial","no sources set up - need sources for initial iteration(s)")
       if(use_exact_nu) then
          n_last_photons_sources = 0
       else
          if(n_last_photons > 0) call error("setup_initial","no sources set up - need sources for last iteration")
       end if
       if(use_raytracing) n_raytracing_photons_sources = 0
    else
       if(use_exact_nu) then
          call mp_read_keyword(input_handle, '/', 'n_last_photons_sources', n_last_photons_sources)
       else
          call mp_read_keyword(input_handle, '/', 'n_last_photons', n_last_photons)
       end if
       if(use_raytracing) call mp_read_keyword(input_handle, '/', 'n_ray_photons_sources', n_raytracing_photons_sources)
    end if

    ! OUTPUT

    g_output = mp_open_group(input_handle, '/Output')

    call mp_read_keyword(g_output, '.', 'output_density', output_density)

    if(trim(output_density).ne.'all' &
         & .and.trim(output_density).ne.'last' &
         & .and.trim(output_density).ne.'none') &
         & call error("setup_initial","output_density should be one of all/last/none")

    call mp_read_keyword(g_output, '.', 'output_density_diff', output_density_diff)

    if(trim(output_density_diff).ne.'all' &
         & .and.trim(output_density_diff).ne.'last' &
         & .and.trim(output_density_diff).ne.'none') &
         & call error("setup_initial","output_density_diff should be one of all/last/none")

    call mp_read_keyword(g_output, '.', 'output_specific_energy', output_specific_energy)

    if(trim(output_specific_energy).ne.'all' &
         & .and.trim(output_specific_energy).ne.'last' &
         & .and.trim(output_specific_energy).ne.'none') &
         & call error("setup_initial","output_specific_energy should be one of all/last/none")

    call mp_read_keyword(g_output, '.', 'output_n_photons', output_n_photons)

    if(trim(output_n_photons).ne.'all' &
         & .and.trim(output_n_photons).ne.'last' &
         & .and.trim(output_n_photons).ne.'none') &
         & call error("setup_initial","output_n_photons should be one of all/last/none")

    call mp_close_group(g_output)

    ! TEMPERATURE CONVERGENCE
    if(n_initial_iter > 0) then
       call mp_read_keyword(input_handle, '/', 'check_convergence', check_convergence)
       if(check_convergence) then
          call mp_read_keyword(input_handle, '/', 'convergence_absolute', convergence_absolute)
          call mp_read_keyword(input_handle, '/', 'convergence_relative', convergence_relative)
          call mp_read_keyword(input_handle, '/', 'convergence_percentile', convergence_percentile)
       end if
    end if

    ! In version 1 dust files there was a bug that caused the Rosseland mean
    ! opacity to be mis-computed (it was in fact computing the Planck inverse
    ! opacity). However, this only affects models that use the PDA, so we only
    ! need to raise an error for these.
    if(use_pda .and. n_dust > 0) then
       do idust=1,n_dust
          if(d(idust)%version == 1) then
             call error("setup_initial", "version 1 dust files can no longer be used when PDA is computed due to a bug - to fix this, re-generate the dust file using the latest version of Hyperion")
          end if
       end do
    end if

  end subroutine setup_initial

  subroutine setup_final_iteration(input_handle)

    implicit none

    integer(hid_t),intent(in) :: input_handle
    integer :: n_peeled
    integer(hid_t) :: g_binned, g_peeled
    character(len=255),allocatable :: group_names(:)

    ! Read configuration for binned images
    g_binned = mp_open_group(input_handle, '/Output/Binned')
    call mp_list_groups(g_binned, '.', group_names)

    if(size(group_names)==0) then
       make_binned_images = .false.
    else if(size(group_names)==1) then
       make_binned_images = .true.
    else
       call error("setup_final_iteration","can't have more than one binned image group")
    end if

    if(make_binned_images) then
       if(use_exact_nu) call error("setup_final_iteration","can't use binned images in exact wavelength mode")
       if(forced_first_interaction) call error("setup_final_iteration", "can't use binned images with forced first interaction")
       call binned_images_setup(g_binned, group_names(1))
    end if

    ! Read configuration for peeloff images
    g_peeled = mp_open_group(input_handle, '/Output/Peeled')
    call mp_list_groups(g_peeled, '.', group_names)
    n_peeled = size(group_names)
    make_peeled_images = n_peeled > 0

    if(make_peeled_images) then
       if(allocated(frequencies)) then
          call peeled_images_setup(g_peeled, group_names, use_raytracing, use_exact_nu, frequencies)
       else
          call peeled_images_setup(g_peeled, group_names, use_raytracing, use_exact_nu)
       end if
    end if

  end subroutine setup_final_iteration

end module setup
