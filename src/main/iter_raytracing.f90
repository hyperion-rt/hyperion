module iteration_raytracing

  use core_lib
  use type_photon
  use sources
  use type_source
  use mpi_core
  use mpi_routines
  use peeled_images
  use dust_main
  use grid_physics

  use grid_geometry

  implicit none
  save

  private
  public :: do_raytracing

contains

  subroutine do_raytracing(n_photons_sources,n_photons_thermal,n_photons_chunk)

    implicit none

    ! Number of photons to run, and size of chunk to use
    integer(idp),intent(in) :: n_photons_sources, n_photons_thermal, n_photons_chunk

    ! Number of photons to run in chunk, and number of photons emitted so far
    integer(idp) :: n_photons, n_photons_curr

    ! Photon object and variable to loop over photons
    integer :: ip
    type(photon) :: p

    n_photons_curr = 0

    ! Tell multi-process routines that this is the start of an iteration
    call mp_reset_first()

    call precompute_jnu_var()

    call mp_join()

    if(main_process()) then
       write(*,*)
       write(*,'("   # Photons    CPU time (sec)    Photons/sec  ")')
       write(*,'(" ----------------------------------------------")')
    end if

    ! Start loop over chunks of photons
    do

       ! Find out how many photons to run
       call mp_n_photons(n_photons_sources, n_photons_curr, n_photons_chunk, n_photons)

       if(n_photons==0) exit

       ! Compute all photons in chunk
       do ip=1,n_photons

          ! Emit photon from original sources
          call emit(p)
          p%energy = p%energy * energy_total / dble(n_photons_sources)
          call peeloff_photon(p, polychromatic=.true.)

       end do

    end do

    call mp_join()

    if(n_dust==0._dp) return

    call mp_reset_first()
    
    if(main_process()) then
       write(*,*)
       write(*,'("   # Photons    CPU time (sec)    Photons/sec  ")')
       write(*,'(" ----------------------------------------------")')
    end if
    
    call mp_join()

    n_photons_curr = 0

    ! Start loop over chunks of photons
    do

       ! Find out how many photons to run
       call mp_n_photons(n_photons_thermal, n_photons_curr, n_photons_chunk, n_photons)

       if(n_photons==0) exit

       ! Compute all photons in chunk
       do ip=1,n_photons

          p = emit_from_grid()

          if(p%energy > 0._dp) then

             p%energy = p%energy * energy_abs_tot(p%dust_id) / dble(n_photons_thermal) * dble(n_dust)

             call peeloff_photon(p, polychromatic=.true.)

          end if

       end do


    end do

    call mp_join()

    ! Need to worry about adding images up with MPI, are we adding the flux up twice?

  end subroutine do_raytracing

end module iteration_raytracing
