module surface_interact

  use type_photon, only : photon
  use type_surface, only : surface
  use surface_collection, only : surfaces
  use sources, only : emit

  implicit none
  save

  private
  public :: interact_with_surface

contains

  subroutine interact_with_surface(p, inu)

    type(photon),intent(inout) :: p

    type(surface), pointer :: srf

    integer,intent(in),optional :: inu

    srf => surfaces(p%surface_id)

    if(srf%source_id > 0) then

       ! Surface is associated with source, re-emit from source

       ! The parentheses are required in the following expression to
       ! force the evaluation of the option (otherwise it gets reset
       ! because p has intent(out) from emit)

       call emit(p, reemit=.true., reemit_id=srf%source_id, reemit_energy=(p%energy), inu=inu)

    else

       ! Scatter

    end if

  end subroutine interact_with_surface

end module surface_interact
