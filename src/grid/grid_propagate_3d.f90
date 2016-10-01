module grid_propagate

  use core_lib
  use type_photon, only : photon
  use type_grid_cell
  use dust_main, only : n_dust
  use grid_geometry, only : escaped, find_wall, in_correct_cell, next_cell, opposite_wall
  use grid_physics, only : specific_energy_sum, density, n_photons, last_photon_id
  use sources
  use counters
  use settings, only : frac_check => propagation_check_frequency

  implicit none
  save

  private
  public :: grid_propagate_debug
  public :: grid_integrate
  public :: grid_integrate_noenergy
  public :: grid_escape_tau
  public :: grid_escape_column_density

  integer :: id

  logical :: debug = .false.

contains

  subroutine grid_propagate_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_propagate_debug

  subroutine grid_integrate(p,tau_required,tau_achieved)

    implicit none

    type(photon),intent(inout) :: p

    real(dp),intent(in)  :: tau_required
    real(dp),intent(out) :: tau_achieved
    real(dp) :: tau_needed,tau_cell,chi_rho_total
    ! the total optical depth to travel, the actual optical depth traveled,
    ! the optical depth remaining to be traveled, and the optical depth in the current cell

    real(dp)  :: tmin, tact
    ! distance to the x, y, and z walls, as well as the closest wall

    type(wall_id) :: id_min

    logical :: radial

    real(dp) :: t_achieved, t_source

    real(dp) :: xi

    integer :: source_id

    radial = (p%r .dot. p%v) > 0.

    if(debug) write(*,'(" [debug] start grid_integrate")')

    ! Set tau achieved to zero

    tau_achieved = 0._dp

    ! Find what cell we are in if we don't know

    if(.not.p%in_cell) call error("grid_integrate", "photon has not been placed in a cell")

    if(allocated(n_photons)) then
       if(last_photon_id(p%icell%ic).ne.p%id) then
          n_photons(p%icell%ic) = n_photons(p%icell%ic) + 1
          last_photon_id(p%icell%ic) = p%id
       end if
    end if

    if(tau_required==0._dp) return

    ! Check what the distance to the nearest source is

    call find_nearest_source(p%r, p%v, t_source, source_id)
    t_achieved = 0._dp

    ! Start the big integral loop

    do

       call random(xi)
       if(xi < frac_check) then
          if(.not.in_correct_cell(p)) then
             call warn("grid_integrate","not in correct cell - killing")
             killed_photons_geo = killed_photons_geo + 1
             p%killed = .true.
             return
          end if
       end if
       ! Find how much still needs to be traveled

       tau_needed = tau_required - tau_achieved

       call find_wall(p,radial,tmin,id_min)

       if(id_min == no_wall) then
          call warn("grid_integrate","cannot find next wall - killing")
          killed_photons_geo = killed_photons_geo + 1
          p%killed = .true.
          return
       end if

       chi_rho_total = 0._dp
       do id=1,n_dust
          chi_rho_total = chi_rho_total + p%current_chi(id) * density(p%icell%ic, id)
       end do
       tau_cell = chi_rho_total * tmin

       if(tau_cell < tau_needed) then

          ! Just move to next cell

          t_achieved = t_achieved + tmin

          if(t_achieved > t_source) then
             p%reabsorbed = .true.
             p%reabsorbed_id = source_id
             return
          end if

          p%r = p%r + tmin * p%v
          tau_achieved = tau_achieved + tau_cell

          ! NEED INDIVIDUAL ALPHA HERE

          do id=1,n_dust
             if(density(p%icell%ic, id) > 0._dp) then
                specific_energy_sum(p%icell%ic, id) = &
                     & specific_energy_sum(p%icell%ic, id) + tmin * p%current_kappa(id) * p%energy
             end if
          end do

          p%on_wall = .true.
          p%icell = next_cell(p%icell, id_min, intersection=p%r)
          p%on_wall_id = opposite_wall(id_min)

          if(debug) write(*,'(" [debug] cross wall")')

          if(escaped(p)) then
             if(debug) write(*,'(" [debug] escape grid")')
             return
          end if

          if(allocated(n_photons)) then
             if(last_photon_id(p%icell%ic).ne.p%id) then
                n_photons(p%icell%ic) = n_photons(p%icell%ic) + 1
                last_photon_id(p%icell%ic) = p%id
             end if
          end if

       else

          ! Move to location in current cell

          tact = tmin * (tau_needed / tau_cell)

          t_achieved = t_achieved + tact

          if(t_achieved > t_source) then
             p%reabsorbed = .true.
             p%reabsorbed_id = source_id
             return
          end if

          p%r = p%r + tact * p%v

          tau_achieved = tau_achieved + tau_needed

          p%on_wall    = .false.
          p%on_wall_id = no_wall

          do id=1,n_dust
             if(density(p%icell%ic, id) > 0._dp) then
                specific_energy_sum(p%icell%ic, id) = &
                     & specific_energy_sum(p%icell%ic, id) &
                     & + tact * p%current_kappa(id) * p%energy
             end if
          end do

          if(debug) write(*,'(" [debug] end grid_integrate")')
          return

       end if

    end do

  end subroutine grid_integrate


  subroutine grid_integrate_noenergy(p,tau_required,tau_achieved)

    implicit none

    type(photon),intent(inout) :: p

    real(dp),intent(in)  :: tau_required
    real(dp),intent(out) :: tau_achieved
    real(dp) :: tau_needed,tau_cell,chi_rho_total
    ! the total optical depth to travel, the actual optical depth traveled,
    ! the optical depth remaining to be traveled, and the optical depth in the current cell

    real(dp)  :: tmin, tact
    ! distance to the x, y, and z walls, as well as the closest wall

    type(wall_id) :: id_min

    logical :: radial

    real(dp) :: t_achieved, t_source

    real(dp) :: xi

    integer :: source_id

    radial = (p%r .dot. p%v) > 0.

    if(debug) write(*,'(" [debug] start grid_integrate_noenergy")')

    ! Set tau achieved to zero

    tau_achieved = 0._dp

    ! Find what cell we are in if we don't know

    if(.not.p%in_cell) call error("grid_integrate_noenergy", "photon has not been placed in a cell")

    if(tau_required==0._dp) return

    ! Check what the distance to the nearest source is

    call find_nearest_source(p%r, p%v, t_source, source_id)
    t_achieved = 0._dp

    ! Start the big integral loop

    do

       call random(xi)
       if(xi < frac_check) then
          if(.not.in_correct_cell(p)) then
             call warn("grid_integrate_noenergy","not in correct cell - killing")
             killed_photons_geo = killed_photons_geo + 1
             p%killed = .true.
             return
          end if
       end if

       ! Find how much still needs to be traveled

       tau_needed = tau_required - tau_achieved

       call find_wall(p,radial,tmin,id_min)

       if(id_min == no_wall) then
          call warn("grid_integrate_noenergy","cannot find next wall - killing")
          killed_photons_geo = killed_photons_geo + 1
          p%killed = .true.
          return
       end if

       chi_rho_total = 0._dp
       do id=1,n_dust
          chi_rho_total = chi_rho_total + p%current_chi(id) * density(p%icell%ic, id)
       end do
       tau_cell = chi_rho_total * tmin

       if(tau_cell < tau_needed) then

          ! Just move to next cell

          t_achieved = t_achieved + tmin

          if(t_achieved > t_source) then
             p%reabsorbed = .true.
             p%reabsorbed_id = source_id
             return
          end if

          p%r = p%r + tmin * p%v
          tau_achieved = tau_achieved + tau_cell

          p%on_wall = .true.
          p%icell = next_cell(p%icell, id_min, intersection=p%r)
          p%on_wall_id = opposite_wall(id_min)

          if(debug) write(*,'(" [debug] cross wall")')

          if(escaped(p)) then
             if(debug) write(*,'(" [debug] escape grid")')
             return
          end if

       else

          ! Move to location in current cell

          tact = tmin * (tau_needed / tau_cell)

          t_achieved = t_achieved + tact

          if(t_achieved > t_source) then
             p%reabsorbed = .true.
             p%reabsorbed_id = source_id
             return
          end if

          p%r = p%r + tact * p%v

          tau_achieved = tau_achieved + tau_needed

          p%on_wall    = .false.
          p%on_wall_id = no_wall

          if(debug) write(*,'(" [debug] end grid_integrate")')
          return

       end if

    end do

  end subroutine grid_integrate_noenergy

  subroutine grid_escape_tau(p_orig,tmax,tau,killed)

    implicit none

    ! This should work with recent compilers but not widespread yet
    ! type(photon),value :: p
    type(photon),intent(in) :: p_orig
    type(photon) :: p

    real(dp),intent(in) :: tmax
    real(dp),intent(out) :: tau
    logical,intent(out) :: killed

    real(dp)  :: tmin, t_source, t_current
    type(wall_id) :: id_min
    logical :: radial, finished
    real(dp) :: xi
    integer :: source_id

    p = p_orig

    killed = .false.

    radial = (p%r .dot. p%v) > 0.

    if(debug) write(*,'(" [debug] start grid_escape_tau")')

    ! Find what cell we are in if we don't know

    if(.not.p%in_cell) call error("grid_escape_tau", "photon has not been placed in a cell")

    ! Check what the distance to the nearest source is

    call find_nearest_source(p%r, p%v, t_source, source_id)

    ! If we are in peeloff mode, there's no point in going on
    if(t_source < tmax) then
       killed = .true.
       return
    end if

    ! Start the big integral loop

    tau = 0._dp
    t_current = 0._dp
    finished = .false.

    do

       call random(xi)
       if(xi < frac_check) then
          if(.not.in_correct_cell(p)) then
             call warn("grid_escape_tau","not in correct cell - killing")
             killed_photons_geo = killed_photons_geo + 1
             killed = .true.
             return
          end if
       end if

       call find_wall(p,radial,tmin,id_min)

       if(id_min == no_wall) then
          call warn("grid_escape_tau","cannot find next wall - killing")
          killed_photons_geo = killed_photons_geo + 1
          killed = .true.
          return
       end if

       if(t_current+tmin > tmax) then
          tmin = tmax - t_current
          finished = .true.
       end if

       p%r = p%r + tmin * p%v
       t_current = t_current + tmin

       do id=1,n_dust
          tau = tau + p%current_chi(id) * density(p%icell%ic, id) * tmin
       end do

       if(finished) return

       p%on_wall = .true.
       p%icell = next_cell(p%icell, id_min, intersection=p%r)
       p%on_wall_id = opposite_wall(id_min)

       if(debug) write(*,'(" [debug] cross wall")')

       if(escaped(p)) then
          if(debug) write(*,'(" [debug] escape grid")')
          return
       end if

    end do

  end subroutine grid_escape_tau

  subroutine grid_escape_column_density(p_orig,tmax,column_density,killed)

    implicit none

    ! This should work with recent compilers but not widespread yet
    ! type(photon),value :: p
    type(photon),intent(in) :: p_orig
    type(photon) :: p

    real(dp),intent(in) :: tmax
    real(dp),intent(out) :: column_density(:)
    logical,intent(out) :: killed

    real(dp)  :: tmin, t_source, t_current
    type(wall_id) :: id_min
    logical :: radial, finished
    real(dp) :: xi
    integer :: source_id

    p = p_orig

    killed = .false.

    radial = (p%r .dot. p%v) > 0.

    if(debug) write(*,'(" [debug] start grid_escape_column_density")')

    ! Find what cell we are in if we don't know

    if(.not.p%in_cell) call error("grid_escape_column_density", "photon has not been placed in a cell")

    ! Check what the distance to the nearest source is

    call find_nearest_source(p%r, p%v, t_source, source_id)

    ! If we are in peeloff mode, there's no point in going on
    if(t_source < tmax) then
       killed = .true.
       return
    end if

    ! Start the big integral loop

    column_density = 0._dp
    t_current = 0._dp
    finished = .false.

    do

       call random(xi)
       if(xi < frac_check) then
          if(.not.in_correct_cell(p)) then
             call warn("grid_escape_column_density","not in correct cell - killing")
             killed_photons_geo = killed_photons_geo + 1
             killed = .true.
             return
          end if
       end if

       call find_wall(p,radial,tmin,id_min)

       if(id_min == no_wall) then
          call warn("grid_escape_column_density","cannot find next wall - killing")
          killed_photons_geo = killed_photons_geo + 1
          killed = .true.
          return
       end if

       if(t_current+tmin > tmax) then
          tmin = tmax - t_current
          finished = .true.
       end if

       p%r = p%r + tmin * p%v
       t_current = t_current + tmin
       column_density = column_density + density(p%icell%ic, :) * tmin

       if(finished) return

       p%on_wall = .true.
       p%icell = next_cell(p%icell, id_min, intersection=p%r)
       p%on_wall_id = opposite_wall(id_min)

       if(debug) write(*,'(" [debug] cross wall")')

       if(escaped(p)) then
          if(debug) write(*,'(" [debug] escape grid")')
          return
       end if

    end do

  end subroutine grid_escape_column_density

end module grid_propagate
