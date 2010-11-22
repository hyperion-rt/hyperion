module grid_pda_geometry

  use core_lib
  use type_grid_cell
  use grid_geometry
  implicit none
  save

  private
  public :: geometrical_factor
  public :: check_allowed_pda

contains

  real(dp) function geometrical_factor(wall, cell)
    implicit none
    integer, intent(in) :: wall
    type(grid_cell), intent(in) :: cell
    geometrical_factor = 1._dp
  end function geometrical_factor

  subroutine check_allowed_pda(do_pda)
    implicit none
    logical,intent(inout) :: do_pda(:)
    do_pda = .false.
  end subroutine check_allowed_pda

end module grid_pda_geometry
