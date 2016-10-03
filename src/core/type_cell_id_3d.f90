module type_grid_cell

  use type_grid

  implicit none
  private

  public :: cell_id
  public :: grid_cell
  type grid_cell
     integer :: i1, i2, i3, ic
  end type grid_cell

  type(grid_cell),parameter,public :: invalid_cell = grid_cell(-1,-1,-1,-1)

  ! Note: if one of the module procedures for .eq. is called equal, this causes
  !       issues with ifort 16 due a bug (because lib_version also defines
  !       equal). Therefore, we call the functions equal_grid_cell and
  !       equall_wall_id instead.

  public :: operator(.eq.)
  interface operator(.eq.)
     module procedure equal_grid_cell
     module procedure equal_wall_id
  end interface operator(.eq.)

  public :: operator(+)
  interface operator(+)
     module procedure add_wall
  end interface operator(+)

  public :: new_grid_cell
  interface new_grid_cell
     module procedure new_grid_cell_1d
     module procedure new_grid_cell_3d
  end interface new_grid_cell

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

  type(wall_id) function add_wall(a,b) result(c)
    implicit none
    type(wall_id), intent(in) :: a,b
    c%w1 = a%w1 + b%w1
    c%w2 = a%w2 + b%w2
    c%w3 = a%w3 + b%w3
  end function add_wall

  logical function equal_grid_cell(a,b)
    implicit none
    type(grid_cell), intent(in) :: a,b
    equal_grid_cell = a%i1 == b%i1 .and. a%i2 == b%i2 .and. a%i3 == b%i3
  end function equal_grid_cell

  type(grid_cell) function new_grid_cell_3d(i1, i2, i3, geo) result(cell)
    implicit none
    integer,intent(in) :: i1, i2, i3
    type(grid_geometry_desc),intent(in) :: geo
    cell%i1 = i1
    cell%i2 = i2
    cell%i3 = i3
    cell%ic = cell_id(i1, i2, i3, geo)
  end function new_grid_cell_3d

  type(grid_cell) function new_grid_cell_1d(ic, geo) result(cell)
    implicit none
    integer,intent(in) :: ic
    type(grid_geometry_desc),intent(in) :: geo
    cell%ic = ic
    cell%i3 = ic / (geo%n1 * geo%n2)
    cell%i2 = (ic - cell%i3 * geo%n1 * geo%n2) / geo%n1
    cell%i1 = ic - cell%i3 * geo%n1 * geo%n2 - cell%i2 * geo%n1
    cell%i2 = cell%i2 + 1
    cell%i3 = cell%i3 + 1
    if(cell%i1==0) then
       cell%i1=geo%n1
       cell%i2=cell%i2 - 1
    end if
    if(cell%i2==0) then
       cell%i2=geo%n2
       cell%i3=cell%i3 - 1
    end if
  end function new_grid_cell_1d

  integer function cell_id(i1, i2, i3, geo)
    implicit none
    integer,intent(in) :: i1, i2, i3
    type(grid_geometry_desc),intent(in) :: geo
    cell_id = (i3-1)*geo%n1*geo%n2 + (i2-1)*geo%n1 + i1
  end function cell_id

end module type_grid_cell
