module dust_interact

  use core_lib

  use type_dust
  use type_photon

  use dust_main
  use grid_physics

  use lorentz, only : doppler_shift

  use type_vector3d

  implicit none
  save

  private
  public :: interact
  public :: interact_peeloff

contains


  subroutine interact(p)

    implicit none

    type(photon),intent(inout) :: p
    integer :: id
    real(dp) :: xi,albedo
    real(dp) :: energy_scaling

    ! given the density and energy of each dust type, process a photon
    ! this means finding out whether to scatter or aborb the photon, and to do
    ! whatever needs doing

    id = select_dust_chi_rho(p)

    albedo = p%current_albedo(id)

    ! Save last direction before scattering/absorbing
    p%a_prev = p%a
    p%v_prev = p%v
    p%s_prev = p%s

    ! Transform frequency to frame of reference of dust
    if(moving) then
       p%nu0 = doppler_shift(p%nu, p%a, vector3d_dp(0._dp, 0._dp, 0._dp)-velocity(p%icell%ic,id))
    else
       p%nu0 = p%nu
    end if

    ! Decide whether to absorb or scatter
    call random(xi)
    if(xi > albedo) then
       call dust_emit(d(id),jnu_var_id(p%icell%ic,id),jnu_var_frac(p%icell%ic,id),p%nu0,p%a,p%s,energy_scaling)
       p%energy = p%energy * energy_scaling
       call update_optconsts(p)
       p%scattered=.false.
       p%reprocessed=.true.
       p%last_isotropic = .true.
       p%dust_id = id
       p%last = 'de'
    else
       call dust_scatter(d(id),p%nu0,p%a,p%s)
       p%scattered=.true.
       p%last_isotropic = .false.
       p%dust_id = id
       p%last = 'ds'
       p%n_scat = p%n_scat + 1
    end if

    ! Lorentz shift into absolute frame of reference
    if(moving) then
       p%nu = doppler_shift(p%nu0, p%a, velocity(p%icell%ic,id))
       p%last_isotropic = .false.  ! otherwise peeloff doesn't get called
    else
       p%nu = p%nu0
    end if

    call angle3d_to_vector3d(p%a,p%v)

  end subroutine interact

  subroutine interact_peeloff(p,a_req)
    implicit none
    type(photon),intent(inout) :: p
    type(angle3d_dp),intent(in)    :: a_req

    select case(p%last)
    case('ds')
       call dust_scatter_peeloff(d(p%dust_id),p%nu0,p%a,p%s,a_req)
    case('de')
       call dust_emit_peeloff(d(p%dust_id),p%nu0,p%a,p%s,a_req)
    case default
       call error("interact_peeloff","unexpected p%last flag: "//p%last)
    end select
    call angle3d_to_vector3d(p%a,p%v)

    ! Lorentz shift here too
    if(moving) then
       p%nu = doppler_shift(p%nu0, p%a, velocity(p%icell%ic,p%dust_id))
    else
       p%nu = p%nu0
    end if

  end subroutine interact_peeloff

end module dust_interact
