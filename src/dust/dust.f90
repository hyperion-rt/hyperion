module dust_main

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use type_dust
  use type_photon
  use settings

  implicit none
  save

  private
  public :: setup_dust
  public :: prepare_photon
  public :: update_optconsts
  public :: chi_planck
  public :: kappa_planck
  public :: chi_inv_planck
  public :: kappa_inv_planck
  public :: chi_rosseland
  public :: kappa_rosseland

  integer,public :: n_dust = 0
  type(dust),allocatable,public :: d(:)

contains

  subroutine setup_dust(group)

    implicit none
    integer(hid_t),intent(in) :: group
    character(len=100),allocatable :: dust_properties(:)
    integer(hid_t) :: g_indiv
    integer :: id

    call mp_list_groups(group, '.', dust_properties)

    n_dust = size(dust_properties)

    allocate(d(n_dust))

    do id=1,n_dust

       if(main_process()) write(*,'(" [dust] reading ",A)') trim(dust_properties(id))

       g_indiv = mp_open_group(group, dust_properties(id))
       call dust_setup(g_indiv,d(id),0._dp)
       call mp_close_group(g_indiv)

    end do

  end subroutine setup_dust

  subroutine prepare_photon(p)
    implicit none
    type(photon),intent(inout) :: p
    allocate(p%current_chi(n_dust))
    allocate(p%current_albedo(n_dust))
    allocate(p%current_kappa(n_dust))
  end subroutine prepare_photon

  subroutine update_optconsts(p)
    implicit none
    type(photon),intent(inout) :: p
    integer :: id
    character(len=1000) :: message
    do id=1,n_dust
       if(p%nu < d(id)%nu(1) .or. p%nu > d(id)%nu(d(id)%n_nu)) then
          write(message, '("photon frequency (",ES10.4," Hz) is outside the range defined for the dust optical properties (",ES10.4," to ",ES10.4," Hz)")') p%nu, d(id)%nu(1), d(id)%nu(d(id)%n_nu)
          call error("update_optconsts", trim(message))
       else
          p%current_chi(id)  = interp1d_loglog(d(id)%nu,d(id)%chi_nu,p%nu)
          p%current_albedo(id) = interp1d_loglog(d(id)%nu,d(id)%albedo_nu,p%nu)
          p%current_kappa(id) = p%current_chi(id) * (1._dp - p%current_albedo(id))
       end if
    end do
  end subroutine update_optconsts

  real(dp) function chi_planck(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    chi_planck = interp1d_loglog(d(id)%specific_energy,d(id)%chi_planck,specific_energy)
  end function chi_planck

  real(dp) function kappa_planck(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    kappa_planck = interp1d_loglog(d(id)%specific_energy,d(id)%kappa_planck,specific_energy)
  end function kappa_planck

  real(dp) function chi_inv_planck(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    chi_inv_planck = interp1d_loglog(d(id)%specific_energy,d(id)%chi_inv_planck,specific_energy)
  end function chi_inv_planck

  real(dp) function kappa_inv_planck(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    kappa_inv_planck = interp1d_loglog(d(id)%specific_energy,d(id)%kappa_inv_planck,specific_energy)
  end function kappa_inv_planck

  real(dp) function chi_rosseland(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    chi_rosseland = interp1d_loglog(d(id)%specific_energy,d(id)%chi_rosseland,specific_energy)
  end function chi_rosseland

  real(dp) function kappa_rosseland(id,specific_energy)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: specific_energy ! specific energy
    kappa_rosseland = interp1d_loglog(d(id)%specific_energy,d(id)%kappa_rosseland,specific_energy)
  end function kappa_rosseland

end module dust_main
