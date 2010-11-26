module iteration_lucy

  use core_lib
  use type_photon
  use sources
  use setup
  use grid_geometry
  use grid_physics
  use mpi_core
  use mpi_routines
  use grid_propagate
  use grid_generic
  use dust_interact
  use dust_main
  use grid_mrw
  use grid_pda
  use settings
  use performance

  implicit none
  save

  private
  public :: do_lucy

contains

  subroutine do_lucy(n_photons_tot,n_photons_chunk)

    implicit none

    ! Number of photons to run, and max size of chunk to use
    integer(idp),intent(in) :: n_photons_tot, n_photons_chunk

    ! Number of photons to run in chunk, and number of photons emitted so far
    integer(idp) :: n_photons, n_photons_curr

    ! Photon object and variable to loop over photons
    integer :: ip, ia
    type(photon) :: p

    ! Number of interactions
    integer :: interactions,mrw_steps

    ! Optical depth to travel
    real(dp) :: tau

    ! REMOVE
    real(dp) :: tau_achieved

    if(n_photons_tot == 0) return

    n_photons_curr = 0

    ! Reset iteration-specific summation variables
    call grid_reset_energy()
    energy_current = 0._dp

    ! Tell multi-process routines that this is the start of an iteration
    call mp_reset_first()

    call precompute_jnu_var()

    if(use_mrw) then
       call update_alpha_rosseland()
       call prepare_mrw()
    end if

    call mp_join()

    if(main_process()) call perf_header()

    ! Start loop over chunks of photons
    do

       ! Find out how many photons to run
       call mp_n_photons(n_photons_tot, n_photons_curr, n_photons_chunk, n_photons)

       if(n_photons==0) exit

       ! Compute all photons in chunk
       do ip=1,n_photons

          ! Emit photon
          call emit(p)

          ! Propagate photon
          do interactions=1, n_inter_max

             if(use_mrw.and.interactions > 1) then
                do mrw_steps=1,n_mrw_max
                   if(tau_rosseland_to_closest_wall(p) > mrw_gamma) then
                      call grid_do_mrw(p)
                   else
                      exit
                   end if
                end do
                if(mrw_steps == n_mrw_max + 1) call error('do_lucy', 'maximum number of MRW steps exceeded')
             end if

             ! Sample a random optical depth and propagate that optical depth
             call random_exp(tau)
             call grid_integrate(p,tau,tau_achieved) 

             if(p%reabsorbed) then

                ! Loop until the photon finally escapes interacting with sources
                do ia=1,n_reabs_max

                   ! The parentheses are required in the following expression to
                   ! force the evaluation of the option (otherwise it gets reset
                   ! because p has intent(out) from emit)
                   call emit(p, reemit=.true., reemit_id=(p%reabsorbed_id))

                   ! Sample optical depth and travel
                   call random_exp(tau)
                   call grid_integrate(p,tau,tau_achieved) 

                   ! If we haven't intersected another source, we can proceed
                   if(.not.p%reabsorbed) exit

                end do

                ! Check that we haven't reached the maximum number of successive reabsorptions
                if(ia == n_reabs_max + 1) call error('do_lucy', 'maximum number of successive re-absorptions exceeded')

             end if

             ! Check whether the photon has escaped the grid or was killed
             if(p%killed.or.escaped(p)) exit

             ! Absorb & re-emit, or scatter
             call interact(p)
             p%killed = kill_on_absorb .and. .not.p%scattered
             if(p%killed) exit

          end do

          if(interactions==n_inter_max+1) then
             call warn("main","photon exceeded maximum number of interactions - killing")
             p%killed = .true.
          end if

       end do

    end do

    call mp_join()

    if(main_process()) call perf_footer()

    ! Collect all summation variables in main process
    call mp_collect()

    ! TEMPORARY
    ! if(main_process()) then
    !   write(filename, '("dtauabs_",I3.3,".fits.gz")') iter
    !   call write_grid(trim(output_dir)//'/'//trim(filename), geo%id, energy_abs)
    ! end if

    if(main_process()) then

       ! Scale the energy so the total is the total requested
       call adjust_energy_grid(energy_total/energy_current)

       ! Compute total energy absorbed for each dust type
       call update_energy_abs_tot()

       ! Compute cell temperatures
       call update_temperature()

       ! Run PDA
       if(use_pda) call solve_pda()

    end if

    ! Broadcast grid arrays out to processes
    call mp_broadcast_temperature()

    ! Sublimate dust if needed
    call sublimate_dust()

  end subroutine do_lucy

end module iteration_lucy
