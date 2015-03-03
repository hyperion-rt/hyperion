module spectrum_gridding

  ! The purpose of this module is to provide a fast way to find the spectrum
  ! of sources binned or interp onto the wavelength grids needed for
  ! images or SEDs. For now, the spectra can be computed once and for all but
  ! in future if we support moving sources then the spectrum will depend on 
  ! the viewing angle. This should be seamless from the outside so this module
  ! provides a way to simply set up the settings for the spectrum
  ! interpolation/binning then to request spectra for a specific viewing angle
  ! and let this module do the rest.

  use sources, only : s
  use core_lib

  implicit none
  save

  private
  public :: set_n_spectral_grids
  public :: set_spectral_grid
  public :: get_spectrum

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
    real(dp),allocatable :: spectrum(:)
  end type cached_spectrum

  type(setting), allocatable :: settings(:)
  type(cached_spectrum), allocatable :: cache(:,:)

contains

  elemental real(dp) function normalized_B_nu(nu,T)
    implicit none
    real(dp),intent(in) :: nu,T
    real(dp),parameter :: a = two * h_cgs / c_cgs / c_cgs / stef_boltz * pi
    real(dp),parameter :: b = h_cgs / k_cgs
    real(dp) :: T4
    T4 = T*T*T*T
    normalized_B_nu = a * nu * nu * nu / ( exp(b*nu/T) - one) / T4
  end function normalized_B_nu

  subroutine set_n_spectral_grids(n_grids)
    implicit none
    integer,intent(in) :: n_grids
    if(allocated(settings)) then
      call error("set_n_spectral_grids", "The number of spectral grids has already been set")
    else
      allocate(settings(n_grids))
      allocate(cache(n_grids, size(s)))
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

  subroutine get_spectrum(source_id, grid_id, spectrum)

    implicit none

    integer,intent(in) :: source_id, grid_id
    real(dp), intent(out) :: spectrum(:)

    ! For now, if the cache is set then we can just use that. In future we
    ! may have to check whether any conditions (e.g. viewing angle) have
    ! changed since last time this was called.

    if(allocated(cache(grid_id, source_id)%spectrum)) then
       spectrum = cache(grid_id, source_id)%spectrum
    else
       if(settings(grid_id)%exact_nu) then
          call get_spectrum_interp(source_id, grid_id, spectrum)
       else
          call get_spectrum_binned(source_id, grid_id, spectrum)
       end if
       allocate(cache(grid_id,source_id)%spectrum(size(spectrum)))
       cache(grid_id,source_id)%spectrum = spectrum
    end if

  end subroutine get_spectrum

  subroutine get_spectrum_interp(source_id, grid_id, spectrum)

    implicit none

    integer,intent(in) :: source_id, grid_id
    real(dp), intent(out) :: spectrum(:)

    select case(s(source_id)%freq_type)
    case(1)
       spectrum = interp1d_loglog(s(source_id)%spectrum%x, s(source_id)%spectrum%pdf, &
            & settings(grid_id)%nu, bounds_error=.false., fill_value=0._dp)
    case(2)
       spectrum = normalized_B_nu(settings(grid_id)%nu,s(source_id)%temperature)
    case default
       print *,"Don't need to set up spectrum for source = ",source_id
    end select

  end subroutine get_spectrum_interp

  subroutine get_spectrum_binned(source_id, grid_id, spectrum)

    implicit none

    integer,intent(in) :: source_id, grid_id
    real(dp), intent(out) :: spectrum(:)
    real(dp),allocatable :: nu(:), fnu(:)
    integer :: inu, n_nu
    real(dp) :: numin, numax
    real(dp) :: log10_nu_min, log10_nu_max

    select case(s(source_id)%freq_type)
    case(1)

       allocate(nu(s(source_id)%spectrum%n))
       allocate(fnu(s(source_id)%spectrum%n))

       nu = s(source_id)%spectrum%x
       fnu = s(source_id)%spectrum%pdf

     case(2)

        log10_nu_min = log10(3.e9_dp)  ! 10 cm (0.05K)
        log10_nu_max = log10(3.e16_dp) ! 10 nm (1/2 million K)

        n_nu = ceiling((log10_nu_max - log10_nu_min) * 100000)
        allocate(nu(n_nu))
        allocate(fnu(n_nu))
        do inu=1,n_nu
           nu(inu) = 10._dp**(real(inu-1,dp)/real(n_nu-1,dp)*(log10_nu_max - log10_nu_min) + log10_nu_min)
        end do

        fnu = normalized_B_nu(nu,s(source_id)%temperature)

    end select

    log10_nu_min = log10(settings(grid_id)%nu_min)
    log10_nu_max = log10(settings(grid_id)%nu_max)

    do inu=1,settings(grid_id)%n_nu

       numin = 10._dp ** (log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu-1, dp) / real(settings(grid_id)%n_nu, dp))
       numax = 10._dp ** (log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu, dp) / real(settings(grid_id)%n_nu, dp))

       spectrum(inu) = integral_loglog(nu, fnu, numin, numax)

    end do

    spectrum = spectrum / integral_loglog(nu, fnu)

  end subroutine get_spectrum_binned

end module spectrum_gridding
