module dust_spectrum_gridding

  ! The purpose of this module is to provide a fast way to find the spectrum
  ! of sources binned or interp onto the wavelength grids needed for
  ! images or SEDs. For now, the spectra can be computed once and for all but
  ! in future if we support moving sources then the spectrum will depend on 
  ! the viewing angle. This should be seamless from the outside so this module
  ! provides a way to simply set up the settings for the spectrum
  ! interpolation/binning then to request spectra for a specific viewing angle
  ! and let this module do the rest.

  use dust_main, only : d
  use core_lib

  implicit none
  save

  private
  public :: set_n_spectral_grids
  public :: set_spectral_grid
  public :: get_dust_spectrum

  interface set_spectral_grid
     module procedure set_spectral_grid_interp
     module procedure set_spectral_grid_binned
  end interface

  type setting

     ! Whether to grid the spectrum by binning (False) or interpolation (True)
     logical :: exact_nu = .false.

     ! For interpolation mode
     real(dp),allocatable :: nu(:)

     ! For binning mode
     integer :: n_nu
     real(dp) :: nu_min, nu_max

     ! Set cache
     real(dp),allocatable :: cache(:,:)

  end type setting
  
  type cached_spectrum
    real(dp),allocatable :: chi
    real(dp),allocatable :: log10_emissivity(:,:)
  end type cached_spectrum

  type(setting), allocatable :: settings(:)
  type(cached_spectrum), target, allocatable :: cache(:,:)

contains

  ! Avoid duplication with source_spectral_gridding module?

  subroutine set_n_spectral_grids(n_grids)
    implicit none
    integer,intent(in) :: n_grids
    if(allocated(settings)) then
      call error("set_n_spectral_grids", "The number of spectral grids has already been set")
    else
      allocate(settings(n_grids))
      allocate(cache(n_grids, size(d)))
    end if
  end subroutine set_n_spectral_grids

  subroutine set_spectral_grid_interp(grid_id, nu)
    implicit none
    integer,intent(in) :: grid_id
    real(dp),intent(in) :: nu(:)
    if(.not.allocated(settings)) call error("set_spectral_grid", "The number of spectral grids has not been set")
    if(grid_id > size(settings)) call error("set_spectral_grid", "The grid_id is larger than the number of settings")
    settings(grid_id)%exact_nu = .true.
    allocate(settings(grid_id)%nu(size(nu)))
    settings(grid_id)%nu = nu
  end subroutine set_spectral_grid_interp

  subroutine set_spectral_grid_binned(grid_id, n_nu, nu_min, nu_max)
    implicit none
    integer,intent(in) :: grid_id
    integer :: n_nu
    real(dp) :: nu_min, nu_max
    if(.not.allocated(settings)) call error("set_spectral_grid", "The number of spectral grids has not been set")
    if(grid_id > size(settings)) call error("set_spectral_grid", "The grid_id is larger than the number of settings")
    settings(grid_id)%exact_nu = .false.
    settings(grid_id)%n_nu = n_nu
    settings(grid_id)%nu_min = nu_min
    settings(grid_id)%nu_max = nu_max
  end subroutine set_spectral_grid_binned

  subroutine get_dust_spectrum(dust_id, grid_id, emiss_var_id, emiss_var_frac, spectrum)

    implicit none

    integer,intent(in) :: dust_id, grid_id
    integer,intent(in) :: emiss_var_id
    real(dp), intent(in) :: emiss_var_frac
    real(dp), intent(out) :: spectrum(:)

    ! For now, if the cache is set then we can just use that. In future we
    ! may have to check whether any conditions (e.g. viewing angle) have
    ! changed since last time this was called, if the dust is moving.

    if(.not.allocated(cache(grid_id, dust_id)%log10_emissivity)) then      
       if(settings(grid_id)%exact_nu) then
          call prepare_spectrum_interp_cache(dust_id, grid_id)
       else
          call prepare_spectrum_binned_cache(dust_id, grid_id)
       end if
    end if

    spectrum = (cache(grid_id, dust_id)%log10_emissivity(:,emiss_var_id + 1) - &
         &      cache(grid_id, dust_id)%log10_emissivity(:,emiss_var_id)) * emiss_var_frac + &
         &      cache(grid_id, dust_id)%log10_emissivity(:,emiss_var_id)
    spectrum = 10._dp**(spectrum)

  end subroutine get_dust_spectrum

  subroutine prepare_spectrum_interp_cache(dust_id, grid_id)

    implicit none

    integer,intent(in) :: dust_id, grid_id
    real(dp), pointer :: e_new(:, :)
    integer :: j

    allocate(cache(grid_id, dust_id)%log10_emissivity(size(settings(grid_id)%nu), size(d(dust_id)%j_nu_var)))

    e_new => cache(grid_id, dust_id)%log10_emissivity

    do j=1,size(d(dust_id)%j_nu_var)
      e_new(:,j) = interp1d_loglog(d(dust_id)%j_nu(j)%x, d(dust_id)%j_nu(j)%pdf, settings(grid_id)%nu, &
                   & bounds_error=.false., fill_value=0._dp)
    end do

  end subroutine prepare_spectrum_interp_cache

  subroutine prepare_spectrum_binned_cache(dust_id, grid_id)

    implicit none

    integer,intent(in) :: dust_id, grid_id
    real(dp),allocatable :: nu(:), fnu(:)
    integer :: inu, n_nu
    real(dp) :: numin, numax
    real(dp) :: log10_nu_min, log10_nu_max
    real(dp), pointer :: e_new(:, :)
    integer :: j

    allocate(cache(grid_id, dust_id)%log10_emissivity(settings(grid_id)%n_nu,size(d(dust_id)%j_nu_var)))

    e_new => cache(grid_id, dust_id)%log10_emissivity

    ! Start off by caching the pre-binned emissivities

    log10_nu_min = log10(settings(grid_id)%nu_min)
    log10_nu_max = log10(settings(grid_id)%nu_max)

    do j=1,size(d(dust_id)%j_nu_var)

       do inu=1,settings(grid_id)%n_nu

          numin = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu-1, dp) / real(settings(grid_id)%n_nu, dp))
          numax = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu, dp) / real(settings(grid_id)%n_nu, dp))

          e_new(inu,j) = integral_loglog(d(dust_id)%j_nu(j)%x, d(dust_id)%j_nu(j)%pdf, numin, numax)

       end do

       e_new(:,j) = e_new(:,j) / integral_loglog(d(dust_id)%j_nu(j)%x, d(dust_id)%j_nu(j)%pdf)

    end do
    
    e_new = log10(e_new)
      
  end subroutine prepare_spectrum_binned_cache

end module dust_spectrum_gridding
