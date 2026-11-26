module mpi_routines

  use core_lib
  use mpi_core
  use performance

  implicit none
  save

  private

  public :: mp_routines_debug
  public :: mp_reset_first
  public :: mp_n_photons
  public :: mp_collect_physical_arrays
  public :: mp_broadcast_specific_energy
  public :: mp_collect_images
  public :: mp_broadcast_convergence
  public :: mp_set_random_seed
  public :: mp_allreduce_min_2d

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

  subroutine mp_routines_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
  end subroutine mp_routines_debug

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

    ! Set a different seed for the random number generator in each process.
    ! The generator uses the absolute value of the seed and gives a degenerate
    ! sequence of random numbers if it is a multiple of 65543 (including
    ! zero), and seeds that differ by a multiple of 65543 give the same
    ! sequence. The seed of each process is therefore moved away from zero as
    ! the rank increases (seed - rank for negative seeds and seed + rank for
    ! positive seeds), so that all processes have different sequences, and an
    ! error is raised if the seed of any process would be degenerate. The
    ! main process (rank 0) uses the seed unchanged.

    implicit none
    integer,intent(in) :: seed
    integer(idp) :: rank_seed

    if(seed == 0) call error("mp_set_random_seed", "seed should not be zero")

    if(seed < 0) then
       rank_seed = int(seed, idp) - int(rank, idp)
    else
       rank_seed = int(seed, idp) + int(rank, idp)
    end if

    if(abs(rank_seed) > int(huge(seed), idp)) then
       call error("mp_set_random_seed", "the seed is too large in absolute value for the number of processes")
    end if

    if(mod(abs(rank_seed), 65543_idp) == 0_idp) then
       call error("mp_set_random_seed", "the seed for one of the processes is a multiple of 65543, which gives a degenerate sequence of random numbers - use a different seed")
    end if

    call set_seed(int(rank_seed))

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


  subroutine mp_allreduce_min_2d(array)
    implicit none
    real(dp),intent(inout) :: array(:,:)
  end subroutine mp_allreduce_min_2d

end module mpi_routines
