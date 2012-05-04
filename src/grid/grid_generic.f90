module grid_generic

  use core_lib
  use mpi_hdf5_io
  use mpi_core

  use grid_io
  use grid_geometry, only : geo
  use grid_physics
  use dust_main
  use settings
  use type_dust

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
  end subroutine grid_reset_energy

  subroutine output_grid(handle, iter, n_iter)

    implicit none

    integer(hid_t),intent(in) :: handle
    integer,intent(in) :: iter, n_iter
    character(len=100) :: group_name
    integer(hid_t) :: group

    if(main_process()) write(*,'(" [output_grid] outputting grid arrays for iteration")')

    write(group_name, '("Iteration ",I5.5)') iter
    group = mp_create_group(handle, group_name)

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
