module surface_interact

  use core_lib, only : vector3d_dp, error, dp, interp1d_loglog
  use type_photon, only : photon
  use type_surface, only : surface, surface_normal, surface_scatter, surface_scatter_peeloff
  use surface_collection, only : surfaces
  use sources, only : emit
  use type_vector3d, only : angle3d_to_vector3d
  use type_angle3d, only : angle3d_dp

  implicit none
  save

  private
  public :: interact_with_surface
  public :: interact_surface_peeloff
contains

  subroutine interact_with_surface(p, inu)

    type(photon),intent(inout) :: p

    type(surface), pointer :: srf

    integer,intent(in),optional :: inu

    type(vector3d_dp) :: n

    real(dp) :: albedo, xi

    srf => surfaces(p%surface_id)

    p%intersected = .false.

    if(srf%source_id > 0) then

       ! Surface is associated with source, re-emit from source

       ! The parentheses are required in the following expression to
       ! force the evaluation of the option (otherwise it gets reset
       ! because p has intent(out) from emit)

       call emit(p, reemit=.true., reemit_id=srf%source_id, reemit_energy=(p%energy), inu=inu)

       p%last = 'sr'

    else

       ! Scatter from source surface

       ! First, interpolate albedo
       albedo = interp1d_loglog(srf%prop%nu, srf%prop%albedo_nu, p%nu)

       ! Sample random value, and either scatter or kill the photon accordingly
       call random(xi)
       if(xi > albedo) then
          p%killed = .true.
       else
          call surface_scatter(srf, p%nu, p%r, p%a, p%s)
          call angle3d_to_vector3d(p%a,p%v)
          p%last = 'su'
       end if

    end if

  end subroutine interact_with_surface

  subroutine interact_surface_peeloff(p, a_req)

    ! Peeloff photon from a surface interaction
    !
    ! Parameters
    ! ----------
    ! p : photon
    !     The photon to peeloff
    ! a_req : angle3d_dp
    !     The requested peeloff angle
    !
    ! Returns
    ! -------
    ! p : photon
    !     The modified photon object

    implicit none

    type(photon),intent(inout) :: p
    type(angle3d_dp),intent(in)    :: a_req

    type(surface), pointer :: srf

    select case(p%last)
    case('su')
       srf => surfaces(p%surface_id)
       call surface_scatter_peeloff(srf,p%nu,p%r, p%a,p%s,a_req)
    case default
       call error("interact_surface_peeloff","unexpected p%last flag: "//p%last)
    end select

    call angle3d_to_vector3d(p%a,p%v)

  end subroutine interact_surface_peeloff

end module surface_interact
