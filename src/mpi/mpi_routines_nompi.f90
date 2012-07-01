module mpi_routines

  use core_lib
  use mpi_core
  use performance

  implicit none
  save

  private

  public :: mp_reset_first
  public :: mp_n_photons
  public :: mp_collect_physical_arrays
  public :: mp_broadcast_specific_energy
  public :: mp_collect_images
  public :: mp_broadcast_convergence
  public :: mp_set_random_seed

  real(dp) :: time_curr
  integer(idp) :: n_photons_chunk

  public :: mp_sync
  interface mp_sync
     module procedure mp_sync_integer4
     module procedure mp_sync_integer8
     module procedure mp_sync_real4
     module procedure mp_sync_real8
  end interface mp_sync

contains

  subroutine mp_reset_first()
    implicit none
    time_curr = 0._dp
    n_photons_chunk = 0
  end subroutine mp_reset_first

  subroutine mp_n_photons(n_photons_tot, n_photons_curr, n_photons_stats, n_photons)
    implicit none

    ! Total number of photons, and size of a chunk
    integer(idp),intent(in) :: n_photons_tot, n_photons_stats

    ! Number of photons requested so far
    integer(idp),intent(inout) :: n_photons_curr

    ! Number of photons to compute this time
    integer(idp),intent(out) :: n_photons

    real(dp),save :: time1 = -1._dp
    real(dp) :: time2

    ! Find the number of photons per chunk
    if(n_photons_chunk == 0) then
       if(n_photons_stats > 0) then
          n_photons_chunk = n_photons_stats
       else
          n_photons_chunk = max(nint(real(n_photons_tot, dp) / 10._dp), 1)
       end if
    end if

    if(time1 < 0._dp) call cpu_time(time1)
    call cpu_time(time2)
    time_curr = time_curr + time2-time1

    if(mod(n_photons_curr,n_photons_chunk)==0.and.n_photons_curr > 0) call perf_numbers(n_photons_curr, time_curr)

    n_photons = min(n_photons_chunk, n_photons_tot - n_photons_curr)
    n_photons_curr = n_photons_curr + n_photons

    time1 = time2

  end subroutine mp_n_photons

  subroutine mp_set_random_seed(seed)
    implicit none
    integer :: seed
    call set_seed(seed + rank)
  end subroutine mp_set_random_seed

  subroutine mp_collect_physical_arrays()
  end subroutine mp_collect_physical_arrays

  subroutine mp_broadcast_specific_energy()
  end subroutine mp_broadcast_specific_energy

  subroutine mp_broadcast_convergence(converged)
    implicit none
    logical,intent(inout) :: converged
  end subroutine mp_broadcast_convergence

  subroutine mp_sync_integer4(value)
    implicit none
    integer,intent(in) :: value
  end subroutine mp_sync_integer4

  subroutine mp_sync_integer8(value)
    implicit none
    integer(idp),intent(in) :: value
  end subroutine mp_sync_integer8

  subroutine mp_sync_real4(value)
    implicit none
    real(sp),intent(in) :: value
  end subroutine mp_sync_real4

  subroutine mp_sync_real8(value)
    implicit none
    real(dp),intent(in) :: value
  end subroutine mp_sync_real8

  subroutine mp_collect_images()
  end subroutine mp_collect_images


end module mpi_routines
