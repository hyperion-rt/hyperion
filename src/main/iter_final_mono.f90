module iteration_final_mono

  use core_lib, only : idp, dp, warn, random_exp, random, error

  use type_photon, only : photon

  use sources, only : emit

  use mpi_core, only : main_process, mp_join
  use mpi_routines, only : mp_reset_first, mp_n_photons

  use peeled_images, only : make_peeled_images, peeloff_photon

  use dust_main, only : n_dust
  use dust_interact, only : interact

  use grid_physics, only : energy_abs_tot, &
       &                   emit_from_grid, &
       &                   precompute_jnu_var

  use grid_monochromatic, only : allocate_monochromatic_grid_pdfs, &
       &                         setup_monochromatic_grid_pdfs, &
       &                         emit_from_monochromatic_grid_pdf, &
       &                         deallocate_monochromatic_grid_pdfs

  use grid_propagate, only : grid_integrate_noenergy, &
       &                     grid_escape_tau

  use grid_geometry, only : escaped

  use settings, only : frequencies, &
       &               n_inter_max, &
       &               n_inter_max_warn, &
       &               n_reabs_max,  &
       &               n_reabs_max_warn,  &
       &               forced_first_interaction, &
       &               forced_first_interaction_algorithm, &
       &               kill_on_scatter, &
       &               monochromatic_energy_threshold

  use performance, only : perf_header, &
       &                  perf_footer

  use counters, only : killed_photons_int

  use forced_interaction, only : forced_interaction_wr99, &
       &                        forced_interaction_baes16, &
       &                         WR99, BAES16

  implicit none
  save

  private
  public :: do_final_mono

contains

  subroutine do_final_mono(n_photons_sources,n_photons_thermal,n_photons_chunk, peeloff_scattering_only)

    implicit none

    ! Number of photons to run, and size of chunk to use
    integer(idp),intent(in) :: n_photons_sources, n_photons_thermal, n_photons_chunk

    ! Whether to only peeloff scattered photons
    logical,intent(in) :: peeloff_scattering_only

    ! Number of photons to run in chunk, and number of photons emitted so far
    integer(idp) :: n_photons, n_photons_curr

    ! Photon object and variable to loop over photons
    integer :: ip
    type(photon) :: p

    integer :: inu
    logical :: empty

    ! Precompute emissivity variable locator for each cell
    call precompute_jnu_var()

    call mp_join()

    if(n_photons_sources > 0) then

       ! Loop over monochromatic frequencies
       do inu=1,size(frequencies)

          if(main_process()) write(*,'(" [mono] computing source photons for nu =",ES11.4," Hz")') frequencies(inu)

          ! Tell multi-process routines that this is the start of an iteration
          call mp_reset_first()

          ! Print performance header
          if(main_process()) call perf_header()

          call mp_join()

          ! Initialize the number of completed photons
          n_photons_curr = 0

          ! Start loop over chunks of photons
          do

             ! Find out how many photons to run
             call mp_n_photons(n_photons_sources, n_photons_curr, n_photons_chunk, n_photons)

             if(n_photons==0) exit

             ! Compute all photons in chunk
             do ip=1,n_photons

                ! Emit photon from a source
                call emit(p,inu=inu)

                ! Scale the energy by the number of photons
                p%energy = p%energy / dble(n_photons_sources)

                ! Peeloff the photons from the star
                if(make_peeled_images) then
                   if(.not.peeloff_scattering_only) call peeloff_photon(p, polychromatic=.false.)
                end if

                ! Propagate until photon is absorbed again
                call propagate(p)

             end do

          end do

          call mp_join()

          ! Print performance footer
          if(main_process()) call perf_footer()

       end do

    else
       if(main_process()) then
          write(*,*)
          write(*,'("      ---------- Skipping source emission ----------")')
          write(*,*)
       end if
    end if

    call mp_join()

    if(n_photons_thermal > 0) then

       ! Allocate emissivity PDFs
       call allocate_monochromatic_grid_pdfs()

       ! Loop over monochromatic frequencies
       do inu=1,size(frequencies)

          ! Pre-compute emissivity grids for each dust type
          call setup_monochromatic_grid_pdfs(inu, empty)

          if(main_process()) write(*,'(" [mono] computing dust photons for nu =",ES11.4," Hz")') frequencies(inu)

          if(.not.empty) then

             ! Tell multi-process routines that this is the start of an iteration
             call mp_reset_first()

             if(main_process()) call perf_header()

             call mp_join()

             ! Initialize the number of completed photons
             n_photons_curr = 0

             ! Start loop over chunks of photons
             do

                ! Find out how many photons to run
                call mp_n_photons(n_photons_thermal, n_photons_curr, n_photons_chunk, n_photons)

                if(n_photons==0) exit

                ! Compute all photons in chunk
                do ip=1,n_photons

                   p = emit_from_monochromatic_grid_pdf(inu)

                   if(p%energy > 0._dp) then

                      ! Scale energy
                      p%energy = p%energy * energy_abs_tot(p%dust_id) / dble(n_photons_thermal) * dble(n_dust)

                      ! Peeloff the photons from the dust
                      if(make_peeled_images) then
                         if(.not.peeloff_scattering_only) call peeloff_photon(p, polychromatic=.false.)
                      end if

                      ! Propagate until photon is absorbed again
                      call propagate(p)

                   end if

                end do

             end do

             call mp_join()

             if(main_process()) call perf_footer()

          else
             if(main_process()) then
                write(*,*)
                write(*,'("      -------- No emission at this frequency -------")')
                write(*,*)
             end if
          end if


       end do

       ! Deallocate emissivity PDFs
       call deallocate_monochromatic_grid_pdfs()

    else
       if(main_process()) then
          write(*,*)
          write(*,'("      ----------- Skipping dust emission -----------")')
          write(*,*)
       end if
    end if

  end subroutine do_final_mono

  subroutine propagate(p)

    implicit none

    type(photon), intent(inout) :: p
    integer(idp) :: interactions
    real(dp) :: tau_achieved, tau, tau_escape
    type(photon) :: p_tmp
    real(dp) :: xi, weight
    logical :: killed
    integer :: ia

    logical, parameter :: force_scatter = .true.
    real(dp) :: energy_initial

    energy_initial = p%energy

    ! Propagate photon
    do interactions=1, n_inter_max+1

       ! If this is the first interaction and the user requested forced
       ! first interaction, we sample tau and modify the photon energy
       ! using a forced first interaction algorith - otherwise we sample
       ! tau the normal way.

       if(interactions == 1 .and. forced_first_interaction) then
          p_tmp = p
          call grid_escape_tau(p_tmp, huge(1._dp), tau_escape, killed)
          if(tau_escape > 1.e-10_dp .and. .not. killed) then
             select case(forced_first_interaction_algorithm)
             case(WR99)
                call forced_interaction_wr99(tau_escape, tau, weight)
             case(BAES16)
                call forced_interaction_baes16(tau_escape, tau, weight)
             case default
                call error("propagate", "Unknown forced first interaction algorithm")
             end select
             p%energy = p%energy * weight
           else
              call random_exp(tau)
          end if
       else
          call random_exp(tau)
       end if

       call grid_integrate_noenergy(p,tau,tau_achieved)

       if(p%reabsorbed) then

          ! Loop until the photon finally escapes interacting with sources
          do ia=1,n_reabs_max

             ! The parentheses are required in the following expression to
             ! force the evaluation of the option (otherwise it gets reset
             ! because p has intent(out) from emit). Technically speaking
             ! this may not be correct because we should not be re-emitting
             ! with the same frequency.
             call emit(p, reemit=.true., reemit_id=(p%reabsorbed_id),&
                  &    reemit_energy=(p%energy), inu=(p%inu))

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
          if(ia == n_reabs_max + 1) then
             if(n_reabs_max_warn) call warn('do_final_mono', 'maximum number of successive re-absorptions exceeded')
             killed_photons_int = killed_photons_int + 1
             p%killed = .true.
             exit
          end if

       end if

       ! Check whether the photon has escaped the grid or was killed
       if(p%killed.or.escaped(p)) exit

       ! We do the following check here rather than after the loop, because if
       ! we limit for example to one interaction, we need to allow propagation,
       ! interaction, propagation.
       if(interactions==n_inter_max+1) then
          if(n_inter_max_warn) call warn("do_final","photon exceeded maximum number of interactions - killing")
          killed_photons_int = killed_photons_int + 1
          p%killed = .true.
          exit
       end if

       ! In the monochromatic algorithm, we always kill the photon when
       ! absorbed, so we might as well force scatterings to get a better
       ! signal-to-noise there.

       call interact(p, force_scatter=force_scatter)
       p%killed = (kill_on_scatter .and. p%scattered) .or. (force_scatter .and. p%energy < energy_initial * monochromatic_energy_threshold)

       if(p%killed) exit
       if(make_peeled_images) call peeloff_photon(p, polychromatic=.false.)

    end do

  end subroutine propagate

end module iteration_final_mono
