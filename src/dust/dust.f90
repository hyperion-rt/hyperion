module dust_main

  use core_lib

  use type_dust
  use type_photon
  use settings

  implicit none
  save

  private
  public :: setup_dust
  public :: prepare_photon
  public :: specific_energy_abs2temperature
  public :: temperature2specific_energy_abs
  public :: update_optconsts
  public :: chi_planck
  public :: kappa_planck
  public :: chi_rosseland
  public :: kappa_rosseland

  integer,public :: n_dust = 0
  type(dust),allocatable,public :: d(:)
  logical,allocatable,public :: is_lte_dust(:)

contains

  subroutine setup_dust(group)

    implicit none
    integer(hid_t),intent(in) :: group
    character(len=100),allocatable :: dust_types(:)
    integer(hid_t) :: g_indiv
    integer :: id

    if(hdf5_path_exists(group, 'Dust types')) then
        call hdf5_table_read_column_auto(group, 'Dust types', 'name', dust_types)
        n_dust = size(dust_types)
    else
        n_dust = 0
    end if
    
    allocate(d(n_dust))
    allocate(is_lte_dust(n_dust))
    
    do id=1,n_dust

       write(*,'(" [dust] reading ",A)') trim(dust_types(id))

       g_indiv = hdf5_open_group(group, dust_types(id))
       call dust_setup(g_indiv,d(id),0._dp)
       call hdf5_close_group(g_indiv)

       is_lte_dust(id) = d(id)%emiss_var == 'T'

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
    do id=1,n_dust
       if(p%nu < d(id)%nu(1) .or. p%nu > d(id)%nu(d(id)%n_nu)) then
          call warn("update_optconsts","photon frequency out of bounds - setting optical constants to zero")
          p%current_chi(id)  = 0._dp
          p%current_albedo(id) = 0._dp
          p%current_kappa(id) = 0._dp
       else
          p%current_chi(id)  = interp1d_loglog(d(id)%nu,d(id)%chi_nu,p%nu)
          p%current_albedo(id) = interp1d_loglog(d(id)%nu,d(id)%albedo_nu,p%nu)
          p%current_kappa(id) = p%current_chi(id) * (1._dp - p%current_albedo(id))
       end if
    end do
  end subroutine update_optconsts

  real(dp) function chi_planck(id,t)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: t ! temperature
    chi_planck = interp1d_loglog(d(id)%T,d(id)%chi_planck,T)
  end function chi_planck

  real(dp) function kappa_planck(id,t)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: t ! temperature
    kappa_planck = interp1d_loglog(d(id)%T,d(id)%kappa_planck,T)
  end function kappa_planck

  real(dp) function chi_rosseland(id,t)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: t ! temperature
    chi_rosseland = interp1d_loglog(d(id)%T,d(id)%chi_rosseland,T)
  end function chi_rosseland

  real(dp) function kappa_rosseland(id,t)
    implicit none
    integer,intent(in)  :: id ! dust type
    real(dp),intent(in) :: t ! temperature
    kappa_rosseland = interp1d_loglog(d(id)%T,d(id)%kappa_rosseland,T)
  end function kappa_rosseland

  real(dp) function specific_energy_abs2temperature(id,specific_energy_abs) result(t)
    ! Requires T^4 * planck opacity to be monotonically increasing (need to check this)
    implicit none
    integer,intent(in)   :: id  ! dust type
    real(dp),intent(in)  :: specific_energy_abs  ! specific energy absorbed
    if(specific_energy_abs > 0._dp) then
       if(specific_energy_abs < d(id)%energy_abs_per_mass(1)) then
          call warn("specific_energy_abs2temperature","cell temperature below minimum allowed")
          t = d(id)%T(1)
       else if(specific_energy_abs > d(id)%energy_abs_per_mass(d(id)%n_t)) then
          call warn("specific_energy_abs2temperature","cell temperature above maximum allowed")
          t = d(id)%T(d(id)%n_t)
       else        
          t = interp1d_loglog(d(id)%energy_abs_per_mass, d(id)%T, specific_energy_abs)
       end if
       t = max(t, minimum_temperature)
    else
       t = minimum_temperature
    end if
  end function specific_energy_abs2temperature

  real(dp) function temperature2specific_energy_abs(id,t) result(specific_energy_abs)
    implicit none
    integer,intent(in)   :: id  ! dust type
    real(dp),intent(in)  :: t  ! temperature
    specific_energy_abs = interp1d_loglog(d(id)%T, d(id)%energy_abs_per_mass, t)
  end function temperature2specific_energy_abs


end module dust_main
