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

  ! log10 of each dust type's nu / chi_nu / albedo_nu tables, and
  ! whether a dust type shares the previous one's frequency grid
  type log_tables
     real(dp),allocatable :: lnu(:), lchi(:), lalb(:)
  end type log_tables
  type(log_tables),allocatable :: lt(:)
  logical,allocatable :: same_nu_prev(:)

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

  subroutine setup_log_tables()
    ! log10 of the frequency, chi and albedo tables of each dust type.
    ! update_optconsts reads these instead of calling log10 four times per
    ! interpolation; the values are the same log10 results, so the interpolated
    ! chi and albedo are bit-identical to interp1d_loglog.
    implicit none
    integer :: id, i, n
    if(allocated(lt)) return
    allocate(lt(n_dust), same_nu_prev(n_dust))
    do id=1,n_dust
       n = d(id)%n_nu
       allocate(lt(id)%lnu(n), lt(id)%lchi(n), lt(id)%lalb(n))
       do i=1,n
          lt(id)%lnu(i) = log10(d(id)%nu(i))
          lt(id)%lchi(i) = 0._dp
          lt(id)%lalb(i) = 0._dp
          if(d(id)%chi_nu(i) > 0._dp) lt(id)%lchi(i) = log10(d(id)%chi_nu(i))
          if(d(id)%albedo_nu(i) > 0._dp) lt(id)%lalb(i) = log10(d(id)%albedo_nu(i))
       end do
       same_nu_prev(id) = .false.
       if(id > 1) then
          if(d(id)%n_nu == d(id-1)%n_nu) same_nu_prev(id) = all(d(id)%nu == d(id-1)%nu)
       end if
    end do
  end subroutine setup_log_tables

  subroutine update_optconsts(p)
    ! One binary search per distinct frequency grid (dust types with
    ! the same grid share it), precomputed log tables, and the bracketing index
    ! cached on the photon (p%inu_cache, valid while p%nu == p%inu_cache_nu) for
    ! dust_scatter.  Same arithmetic as interp1d_loglog / interp1d_single_loglog.
    implicit none
    type(photon),intent(inout) :: p
    integer :: id, i
    real(dp) :: lnu, frac
    character(len=1000) :: message
    if(.not.allocated(lt)) call setup_log_tables()
    lnu = log10(p%nu)
    i = -1
    do id=1,n_dust
       if(p%nu < d(id)%nu(1) .or. p%nu > d(id)%nu(d(id)%n_nu)) then
          write(message, '("photon frequency (",ES10.4," Hz) is outside the range defined for the dust optical properties (",ES10.4," to ",ES10.4," Hz)")') p%nu, d(id)%nu(1), d(id)%nu(d(id)%n_nu)
          call error("update_optconsts", trim(message))
       else
          if(.not.(id > 1 .and. same_nu_prev(id) .and. i > 0)) i = locate(d(id)%nu, p%nu)
          frac = ( lnu - lt(id)%lnu(i) ) / ( lt(id)%lnu(i+1) - lt(id)%lnu(i) )
          if(d(id)%chi_nu(i)==0. .or. d(id)%chi_nu(i+1)==0.) then
             p%current_chi(id) = 0._dp
          else
             p%current_chi(id) = 10._dp**(lt(id)%lchi(i) + frac * ( lt(id)%lchi(i+1) - lt(id)%lchi(i) ))
          end if
          if(d(id)%albedo_nu(i)==0. .or. d(id)%albedo_nu(i+1)==0.) then
             p%current_albedo(id) = 0._dp
          else
             p%current_albedo(id) = 10._dp**(lt(id)%lalb(i) + frac * ( lt(id)%lalb(i+1) - lt(id)%lalb(i) ))
          end if
          p%current_kappa(id) = p%current_chi(id) * (1._dp - p%current_albedo(id))
          if(id <= n_inu_cache) p%inu_cache(id) = i
       end if
    end do
    p%inu_cache_nu = p%nu
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
