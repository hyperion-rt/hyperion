module dust_interact

  use core_lib

  use type_dust
  use type_photon

  use dust_main
  use grid_physics


  implicit none
  save

  private
  public :: interact
  public :: interact_peeloff

contains


  subroutine interact(p, force_scatter)

    implicit none

    type(photon),intent(inout) :: p
    integer :: id
    real(dp) :: xi,albedo
    real(dp) :: energy_scaling
    logical,optional :: force_scatter

    ! given the density and energy of each dust type, process a photon
    ! this means finding out whether to scatter or aborb the photon, and to do
    ! whatever needs doing

    id = select_dust_chi_rho(p)

    albedo = p%current_albedo(id)

    ! Save last direction before scattering/absorbing
    p%a_prev = p%a
    p%v_prev = p%v
    p%s_prev = p%s

    ! Decide whether to absorb or scatter
    if(present(force_scatter) .and. force_scatter) then
      xi = 0._dp
    else
      call random(xi)
    end if
    
    if(xi > albedo) then
       call dust_emit(d(id),jnu_var_id(p%icell%ic,id),jnu_var_frac(p%icell%ic,id),p%nu,p%a,p%s,energy_scaling)
       p%energy = p%energy * energy_scaling
       call update_optconsts(p)
       p%scattered=.false.
       p%reprocessed=.true.
       p%last_isotropic = .true.
       p%dust_id = id
       p%last = 'de'
    else
       call dust_scatter(d(id),p%nu,p%a,p%s)
       p%scattered=.true.
       p%last_isotropic = .false.
       p%dust_id = id
       p%last = 'ds'
       p%n_scat = p%n_scat + 1
    end if

    call angle3d_to_vector3d(p%a,p%v)

    if(present(force_scatter) .and. force_scatter) then
      p%energy = p%energy * albedo
    end if

  end subroutine interact

  subroutine interact_peeloff(p,a_req)
    implicit none
    type(photon),intent(inout) :: p
    type(angle3d_dp),intent(in)    :: a_req
    select case(p%last)
    case('ds')
       call dust_scatter_peeloff(d(p%dust_id),p%nu,p%a,p%s,a_req)
    case('de')
       call dust_emit_peeloff(d(p%dust_id),p%nu,p%a,p%s,a_req)
    case default
       call error("interact_peeloff","unexpected p%last flag: "//p%last)
    end select
    call angle3d_to_vector3d(p%a,p%v)
  end subroutine interact_peeloff

end module dust_interact
