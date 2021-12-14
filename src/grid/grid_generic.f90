module grid_generic

  use core_lib, only : sp, dp, hid_t, warn, error
  use mpi_hdf5_io, only : mp_create_group
  use mpi_core, only : main_process
  
  use grid_io, only : write_grid_3d, write_grid_4d, write_grid_5d
  use grid_geometry, only : geo
  use grid_physics, only : n_photons, last_photon_id, specific_energy_sum, specific_energy_sum_nu, specific_energy, density, density_original
  use settings, only : output_n_photons, output_specific_energy, output_density, output_density_diff, physics_io_type, compute_isrf

  !DN Crazy Additions
  use dust_main, only: d
  

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
    if(allocated(specific_energy_sum_nu)) specific_energy_sum_nu = 0._dp
  end subroutine grid_reset_energy

  subroutine output_grid(group, iter, n_iter)

    implicit none

    integer(hid_t),intent(in) :: group
    integer,intent(in) :: iter, n_iter
    !DN CRAZY ADDITION
    integer :: n_cells, n_dust, n_isrf_lam
    integer :: i,j,k
    real, dimension(d(1)%n_nu) :: energy_frequency_bins

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


    
    ! DN Crazy Additions
    ! WRITE THE ISRF
    if (compute_isrf) then 
       
       n_cells = size(specific_energy_sum_nu, 1)
       n_dust = size(specific_energy_sum_nu, 2)
       n_isrf_lam = size(specific_energy_sum_nu,3)
       
       do i=1,d(1)%n_nu
          energy_frequency_bins(i) = d(1)%nu(i)
       end do
       
       if(trim(output_specific_energy)=='all' .or. (trim(output_specific_energy)=='last'.and.iter==n_iter)) then
          !if(allocated(energy_frequency_bins)) then
          call write_grid_3d(group, 'ISRF_frequency_bins',energy_frequency_bins, geo)
       else
          call warn("output_grid","energy_frequency bins [ISRF wavelengths] array is not allocated")
          !end if
       end if
       
       if(trim(output_specific_energy)=='all' .or. (trim(output_specific_energy)=='last'.and.iter==n_iter)) then
          if(allocated(specific_energy_sum_nu)) then
             
             select case(physics_io_type)
                
             case(sp)
                call write_grid_5d(group, 'specific_energy_sum_nu', real(specific_energy_sum_nu, sp), geo)
             case(dp)
                call write_grid_5d(group, 'specific_energy_sum_nu', real(specific_energy_sum_nu, dp), geo)
             case default
                call error("output_grid","unexpected value of physics_io_type (should be sp or dp)")
             end select
          else
             call warn("output_grid","specific_energy_sum_nu array is not allocated")
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
