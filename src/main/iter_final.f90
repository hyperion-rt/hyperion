module iteration_final

  use core_lib
  use type_photon
  use sources
  use setup
  use grid_geometry
  use grid_physics
  use mpi_core
  use mpi_routines
  use grid_propagate
  use dust_interact
  use peeled_images
  use binned_images
  use grid_mrw
  use settings
  use performance

  implicit none
  save

  private
  public :: do_final

contains

  subroutine do_final(n_photons_tot,n_photons_chunk,peeloff_scattering_only)

    implicit none

    ! Number of photons to run, and size of chunk to use
    integer(idp),intent(in) :: n_photons_tot, n_photons_chunk

    ! Number of photons to run in chunk, and number of photons emitted so far
    integer(idp) :: n_photons, n_photons_curr

    ! Photon object and variable to loop over photons
    integer :: ip
    type(photon) :: p

    ! Whether to only peeloff scattered photons
    logical,intent(in) :: peeloff_scattering_only

    if(n_photons_tot == 0) then
       if(main_process()) then
          write(*,*)
          write(*,'("      ------------------ Skipping ------------------")')
          write(*,*)
       end if
       return
    end if

    n_photons_curr = 0

    energy_current = 0._dp

    ! Tell multi-process routines that this is the start of an iteration
    call mp_reset_first()

    if(use_mrw) then
       call update_alpha_rosseland()
       call prepare_mrw()
    end if

    call precompute_jnu_var()

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

          ! Peeloff
          if(make_peeled_images) then
             if(.not.peeloff_scattering_only) call peeloff_photon(p, polychromatic=.false.)
          end if

          ! Propagate until photon is absorbed again
          call propagate(p,peeloff_scattering_only)

          ! Bin photon in image (last iteration)
          if(make_binned_images) then
             if(.not.p%killed) call binned_images_bin_photon(p)
          end if

       end do

    end do

    call mp_join()

    if(main_process()) call perf_footer()

    call mp_sync_energy()

    if(make_binned_images) call binned_images_adjust_scale(energy_total/energy_current)
    if(make_peeled_images) call peeled_images_adjust_scale(energy_total/energy_current)

  end subroutine do_final

  subroutine propagate(p,peeloff_scattering_only)

    implicit none

    type(photon), intent(inout) :: p
    logical,intent(in) :: peeloff_scattering_only
    integer(idp) :: interactions
    real(dp) :: tau_achieved, tau, tau_escape
    type(photon) :: p_tmp
    real(dp) :: xi
    logical :: killed
    integer :: mrw_steps
    integer :: ia

    ! Propagate photon
    do interactions=1, n_inter_max

       ! Sample a random optical depth and propagate that optical depth
       call random_exp(tau)

       if(interactions==1) then
          if(forced_first_scattering) then
             p_tmp = p
             call grid_escape_tau(p_tmp, huge(1._dp), tau_escape, killed)
             if(tau_escape > 1.e-10_dp) then
                call random(xi)
                tau = -log(1._dp-xi*(1._dp - exp(-tau_escape)))
                p%energy = p%energy * (1._dp - exp(-tau_escape))
             end if
          end if
       else
          if(use_mrw) then
             do mrw_steps=1,n_mrw_max
                if(tau_rosseland_to_closest_wall(p) > mrw_gamma) then
                   call grid_do_mrw_noenergy(p)
                   if(make_peeled_images) then
                      if(.not.peeloff_scattering_only) call peeloff_photon(p, polychromatic=.false.)
                   end if
                else
                   exit
                end if
             end do
             if(mrw_steps == n_mrw_max+1) call error('do_lucy', 'maximum number of MRW steps exceeded')
          end if
       end if

       call grid_integrate_noenergy(p,tau,tau_achieved) 

       if(p%reabsorbed) then

          ! Loop until the photon finally escapes interacting with sources
          do ia=1,n_reabs_max

             ! The parentheses are required in the following expression to
             ! force the evaluation of the option (otherwise it gets reset
             ! because p has intent(out) from emit)
             call emit(p, reemit=.true., reemit_id=(p%reabsorbed_id))

             ! We now peeloff the photon even if only scattered photons are
             ! wanted because this is a kind of scattering, and will not be
             ! taken into account in the raytracing.
             if(make_peeled_images) call peeloff_photon(p, polychromatic=.false.)

             ! Sample optical depth and travel
             call random_exp(tau)
             call grid_integrate_noenergy(p,tau,tau_achieved) 

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

       ! Peeloff
       if(make_peeled_images) then
          if(p%scattered.or..not.peeloff_scattering_only) call peeloff_photon(p, polychromatic=.false.)
       end if

    end do

    if(interactions==n_inter_max+1) then
       call warn("main","photon exceeded maximum number of interactions - killing")
       p%killed = .true.
    end if

  end subroutine propagate

end module iteration_final
