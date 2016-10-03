module type_grid_cell

  use type_grid

  implicit none
  private

  public :: grid_cell
  type grid_cell
     integer :: ic
  end type grid_cell

  type(grid_cell),parameter,public :: invalid_cell = grid_cell(-1)

  ! Note: if one of the module procedures for .eq. is called equal, this causes
  !       issues with ifort 16 due a bug (because lib_version also defines
  !       equal). Therefore, we call the functions equal_grid_cell and
  !       equall_wall_id instead.

  public :: operator(.eq.)
  interface operator(.eq.)
     module procedure equal_grid_cell
     module procedure equal_wall_id
  end interface operator(.eq.)

  public :: new_grid_cell

  public :: wall_id
  type wall_id
     integer :: w1=0, w2=0, w3=0
  end type wall_id

  type(wall_id), parameter, public :: no_wall = wall_id(0, 0, 0)

contains

  logical function equal_wall_id(a,b)
    implicit none
    type(wall_id), intent(in) :: a,b
    equal_wall_id = a%w1 == b%w1 .and. a%w2 == b%w2 .and. a%w3 == b%w3
  end function equal_wall_id

  logical function equal_grid_cell(a,b)
    implicit none
    type(grid_cell), intent(in) :: a,b
    equal_grid_cell = a%ic == b%ic
  end function equal_grid_cell

  type(grid_cell) function new_grid_cell(ic, geo) result(cell)
    implicit none
    integer,intent(in) :: ic
    type(grid_geometry_desc),intent(in) :: geo
    cell%ic = ic
  end function new_grid_cell

end module type_grid_cell
