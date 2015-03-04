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


end module spectrum_gridding
