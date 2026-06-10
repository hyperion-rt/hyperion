module grid_generic

  use core_lib, only : sp, dp, hid_t, warn, error
  use mpi_hdf5_io, only : mp_create_group, mp_write_array
  use mpi_core, only : main_process
  
  use grid_io, only : write_grid_3d, write_grid_4d, write_grid_5d
  use grid_geometry, only : geo
  use grid_physics, only : n_photons, last_photon_id, specific_energy_sum, specific_energy_sum_spectrum, specific_energy, specific_energy_spectrum, nu_bins, density, density_original
  use settings, only : output_n_photons, output_specific_energy, output_specific_energy_spectrum, output_density, output_density_diff, physics_io_type, compute_specific_energy_spectrum

  implicit none
  save

  private
  public :: grid_reset_energy
  public :: output_grid

contains

  subroutine grid_reset_energy()
    implicit none
    if(allocated(n_photons)) n_photons = 0
    if(allocated(last_photon_id)) last_photon_id = 0
    if(allocated(specific_energy_sum)) specific_energy_sum = 0._dp
    if(allocated(specific_energy_sum_spectrum)) specific_energy_sum_spectrum = 0._dp
  end subroutine grid_reset_energy

  subroutine output_grid(group, iter, n_iter)

    implicit none

    integer(hid_t),intent(in) :: group
    integer,intent(in) :: iter, n_iter

    if(main_process()) write(*,'(" [output_grid] outputting grid arrays for iteration")')

    ! NUMBER OF PHOTONS IN EACH CELL

    if(trim(output_n_photons)=='all' .or. (trim(output_n_photons)=='last'.and.iter==n_iter)) then
       if(allocated(n_photons)) then
          call write_grid_3d(group, 'n_photons', n_photons, geo)
       else
          call warn("output_grid","n_photons array is not allocated")
       end if
    end if

    ! ENERGY/PATH LENGTHS

    if(trim(output_specific_energy)=='all' .or. (trim(output_specific_energy)=='last'.and.iter==n_iter)) then
       if(allocated(specific_energy)) then
          select case(physics_io_type)
          case(sp)
             call write_grid_4d(group, 'specific_energy', real(specific_energy, sp), geo)
          case(dp)
             call write_grid_4d(group, 'specific_energy', real(specific_energy, dp), geo)
          case default
             call error("output_grid","unexpected value of physics_io_type (should be sp or dp)")
          end select
       else
          call warn("output_grid","specific_energy array is not allocated")
       end if
    end if


    
    ! WRITE THE FREQUENCY-RESOLVED SPECIFIC ENERGY
    if (compute_specific_energy_spectrum) then

       if(trim(output_specific_energy_spectrum)=='all' .or. (trim(output_specific_energy_spectrum)=='last'.and.iter==n_iter)) then

          ! The frequencies are a per-frequency quantity, not per-cell, so they
          ! are written as a plain 1-D dataset rather than a grid array.
          call mp_write_array(group, 'specific_energy_spectrum_frequencies', nu_bins)

          if(allocated(specific_energy_spectrum)) then
             select case(physics_io_type)
             case(sp)
                call write_grid_5d(group, 'specific_energy_spectrum', real(specific_energy_spectrum, sp), geo)
             case(dp)
                call write_grid_5d(group, 'specific_energy_spectrum', real(specific_energy_spectrum, dp), geo)
             case default
                call error("output_grid","unexpected value of physics_io_type (should be sp or dp)")
             end select
          else
             call warn("output_grid","specific_energy_spectrum array is not allocated")
          end if

       end if

    endif


    ! DENSITY

    if(trim(output_density)=='all' .or. (trim(output_density)=='last'.and.iter==n_iter)) then
       if(allocated(density)) then
          select case(physics_io_type)
          case(sp)
             call write_grid_4d(group, 'density', real(density, sp), geo)
          case(dp)
             call write_grid_4d(group, 'density', real(density, dp), geo)
          case default
             call error("output_grid","unexpected value of physics_io_type (should be sp or dp)")
          end select
       else
          call warn("output_grid","density array is not allocated")
       end if
    end if


    ! DENSITY DIFFERENCE

    if(trim(output_density_diff)=='all' .or. (trim(output_density_diff)=='last'.and.iter==n_iter)) then
       if(allocated(density).and.allocated(density_original)) then
          select case(physics_io_type)
          case(sp)
             call write_grid_4d(group, 'density_diff', real(density - density_original, sp), geo)
          case(dp)
             call write_grid_4d(group, 'density_diff', real(density - density_original, dp), geo)
          case default
             call error("output_grid","unexpected value of physics_io_type (should be sp or dp)")
          end select
       else
          if(.not.allocated(density)) call warn("output_grid","density array is not allocated")
          if(.not.allocated(density_original)) call warn("output_grid","density_original array is not allocated")
       end if
    end if

  end subroutine output_grid

end module grid_generic
