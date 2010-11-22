module type_grid_cell

  use type_grid

  implicit none
  private

  public :: grid_cell
  type grid_cell
     integer :: ic
  end type grid_cell

  type(grid_cell),parameter,public :: invalid_cell = grid_cell(-1)

  public :: operator(.eq.)
  interface operator(.eq.)
     module procedure equal
  end interface operator(.eq.)

  public :: new_grid_cell

contains

  logical function equal(a,b)
    implicit none
    type(grid_cell), intent(in) :: a,b
    equal = a%ic == b%ic
  end function equal

  type(grid_cell) function new_grid_cell(ic, geo) result(cell)
    implicit none
    integer,intent(in) :: ic
    type(grid_geometry_desc),intent(in) :: geo
    cell%ic = ic
  end function new_grid_cell

end module type_grid_cell
