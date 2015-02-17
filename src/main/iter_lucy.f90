module iteration_lucy

  use core_lib, only : idp, dp, random_exp, random, error, warn

  use type_photon, only : photon

  use sources, only : emit, energy_current, energy_total

  use mpi_core, only : main_process, mp_join
  use mpi_routines, only : mp_reset_first, &
       &                   mp_n_photons, &
       &                   mp_sync, &
       &                   mp_collect_physical_arrays, &
       &                   mp_broadcast_specific_energy

  use peeled_images, only : make_peeled_images, &
       &                    peeloff_photon, &
       &                    peeled_images_adjust_scale
  use binned_images, only : make_binned_images, &
       &                    binned_images_bin_photon, &
       &                    binned_images_adjust_scale

  use dust_interact, only : interact

  use grid_geometry, only : escaped

  use grid_generic, only : grid_reset_energy

  use grid_physics, only : tau_inv_planck_to_closest_wall, &
       &                   precompute_jnu_var, &
       &                   update_alpha_inv_planck, &
       &                   update_energy_abs, &
       &                   sublimate_dust

  use grid_propagate, only : grid_escape_tau, &
       &                     grid_integrate

  use grid_mrw, only : prepare_mrw, &
       &               grid_do_mrw

  use grid_pda, only : solve_pda

  use settings, only : forced_first_scattering, &
       &               kill_on_absorb, kill_on_scatter, &
       &               n_inter_max, &
       &               n_inter_max_warn, &
       &               mrw_gamma, &
       &               n_mrw_max, &
       &               n_mrw_max_warn, &
       &               n_reabs_max, &
       &               n_reabs_max_warn, &
       &               use_mrw, &
       &               use_pda

  use performance, only : perf_header, perf_footer

  use counters, only : killed_photons_int

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

    if(n_photons_tot == 0) then
       if(main_process()) then
          write(*,*)
          write(*,'("      ------------------ Skipping ------------------")')
          write(*,*)
       end if
       return
    end if

    n_photons_curr = 0

    ! Reset iteration-specific summation variables
    call grid_reset_energy()
    energy_current = 0._dp

    ! Tell multi-process routines that this is the start of an iteration
    call mp_reset_first()

    call precompute_jnu_var()

    if(use_mrw) then
       call update_alpha_inv_planck()
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
          do interactions=1, n_inter_max+1

             if(use_mrw.and.interactions > 1) then
                do mrw_steps=1,n_mrw_max
                   if(tau_inv_planck_to_closest_wall(p) > mrw_gamma) then
                      call grid_do_mrw(p)
                   else
                      exit
                   end if
                end do
                if(mrw_steps == n_mrw_max + 1) then
                   if(n_mrw_max_warn) call warn("do_lucy","maximum number of MRW steps exceeded - killing")
                   killed_photons_int = killed_photons_int + 1
                   p%killed = .true.
                   exit
                end if
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
                   call emit(p, reemit=.true., reemit_id=(p%reabsorbed_id), reemit_energy=(p%energy))

                   ! Sample optical depth and travel
                   call random_exp(tau)
                   call grid_integrate(p,tau,tau_achieved)

                   ! If we haven't intersected another source, we can proceed
                   if(.not.p%reabsorbed) exit

                end do

                ! Check that we haven't reached the maximum number of successive reabsorptions
                if(ia == n_reabs_max + 1) then
                   if(n_reabs_max_warn) call warn('do_lucy', 'maximum number of successive re-absorptions exceeded')
                   killed_photons_int = killed_photons_int + 1
                   p%killed = .true.
                   exit
                end if

             end if

             ! Check whether the photon has escaped the grid or was killed
             if(p%killed.or.escaped(p)) exit

             ! We do the following check here rather than after the loop, 
             ! because if we limit for example to one interaction, we need to 
             ! allow propagation, interaction, propagation.
             if(interactions==n_inter_max + 1) then
                if(n_inter_max_warn) call warn("do_lucy","photon exceeded maximum number of interactions - killing")
                killed_photons_int = killed_photons_int + 1
                p%killed = .true.
                exit
             end if

             ! Absorb & re-emit, or scatter
             call interact(p)
             p%killed = (kill_on_scatter .and. p%scattered) .or. (kill_on_absorb .and. .not.p%scattered)
             if(p%killed) exit

          end do

       end do

    end do

    call mp_join()

    if(main_process()) call perf_footer()

    ! Collect all summation arrays in main process
    call mp_collect_physical_arrays()

    ! Sync energy emitted
    call mp_sync(energy_current)

    if(main_process()) then

       ! Scale the energy so the total is the total requested
       call update_energy_abs(energy_total/energy_current)

       ! Run PDA
       if(use_pda) call solve_pda()

    end if

    ! Broadcast specific energy out to processes
    call mp_broadcast_specific_energy()

    ! Sublimate dust if needed
    call sublimate_dust()

  end subroutine do_lucy

end module iteration_lucy
