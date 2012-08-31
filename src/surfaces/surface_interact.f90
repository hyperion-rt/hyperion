module surface_interact

  use core_lib, only : vector3d_dp
  use type_photon, only : photon
  use type_surface, only : surface, surface_normal, surface_scatter
  use surface_collection, only : surfaces
  use sources, only : emit
  use type_vector3d, only : angle3d_to_vector3d

  implicit none
  save

  private
  public :: interact_with_surface

contains

  subroutine interact_with_surface(p, inu)

    type(photon),intent(inout) :: p

    type(surface), pointer :: srf

    integer,intent(in),optional :: inu

    type(vector3d_dp) :: n

    srf => surfaces(p%surface_id)

    if(srf%source_id > 0) then

       ! Surface is associated with source, re-emit from source

       ! The parentheses are required in the following expression to
       ! force the evaluation of the option (otherwise it gets reset
       ! because p has intent(out) from emit)

       call emit(p, reemit=.true., reemit_id=srf%source_id, reemit_energy=(p%energy), inu=inu)

    else

       ! Scatter from source surface

       call surface_scatter(srf, p%nu, p%r, p%a, p%s)
       call angle3d_to_vector3d(p%a,p%v)

       p%intersected = .false.

    end if

  end subroutine interact_with_surface

end module surface_interact
