module mpi_routines

  use core_lib
  use mpi_core
  use performance

  implicit none
  save

contains

  subroutine mp_reset_first()
  end subroutine mp_reset_first

  subroutine mp_n_photons(n_photons_tot, n_photons_curr, n_photons_chunk, n_photons)
    implicit none

    ! Total number of photons, and size of a chunk
    integer(idp),intent(in) :: n_photons_tot, n_photons_chunk

    ! Number of photons requested so far
    integer(idp),intent(inout) :: n_photons_curr

    ! Number of photons to compute this time
    integer(idp),intent(out) :: n_photons

    real(dp) :: time1=-1.
    real(dp) :: time2
    real(dp),save :: time_curr = 0._dp

    if(time1 < 0._dp) call cpu_time(time1)
    call cpu_time(time2)
    time_curr = time_curr + time2-time1

    if(mod(n_photons_curr,n_photons_chunk)==0) call perf_numbers(n_photons_curr, time_curr)

    n_photons = min(n_photons_chunk, n_photons_tot - n_photons_curr)
    n_photons_curr = n_photons_curr + n_photons

    time1 = time2

  end subroutine mp_n_photons

  subroutine mp_set_random_seed()
    implicit none
    call set_seed(-124902+rank)
  end subroutine mp_set_random_seed

  subroutine mp_collect()
  end subroutine mp_collect

  subroutine mp_broadcast_temperature() 
  end subroutine mp_broadcast_temperature

  subroutine mp_sync_energy() 
  end subroutine mp_sync_energy

  subroutine mp_sync_cputime(cputime)
    implicit none
    real(dp),intent(inout) :: cputime
  end subroutine mp_sync_cputime

  subroutine mp_collect_results()
  end subroutine mp_collect_results


end module mpi_routines
