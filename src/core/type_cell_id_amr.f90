module type_grid_cell

  use type_grid

  implicit none
  private

  public :: preset_cell_id

  public :: grid_cell
  type grid_cell
     integer :: ic         ! unique cell ID
     integer :: ilevel     ! level ID
     integer :: igrid       ! grid ID
     integer :: i1, i2, i3 ! position in grid
  end type grid_cell

  type(grid_cell),parameter,public :: invalid_cell = grid_cell(-1, -1, -1, -1, -1, -1)
  type(grid_cell),parameter,public :: outside_cell = grid_cell(-2, -2, -2, -2, -2, -2)

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
  interface new_grid_cell
     module procedure new_grid_cell_1d
     module procedure new_grid_cell_5d
  end interface new_grid_cell

  ! level, grid, and coordinates for each unique ID
  integer,allocatable :: cell_ilevel(:), cell_igrid(:), cell_i1(:), cell_i2(:), cell_i3(:)

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

  subroutine preset_cell_id(geo)

    implicit none

    type(grid_geometry_desc),intent(in), target :: geo
    type(level_desc), pointer :: level
    type(grid_desc), pointer :: grid
    integer :: ic, i1, i2, i3, ilevel, igrid

    allocate(cell_ilevel(geo%n_cells))
    allocate(cell_igrid(geo%n_cells))
    allocate(cell_i1(geo%n_cells))
    allocate(cell_i2(geo%n_cells))
    allocate(cell_i3(geo%n_cells))

    ic = 0
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          do i3=1,grid%n3
             do i2=1,grid%n2
                do i1=1,grid%n1
                   ic = ic + 1
                   cell_ilevel(ic) = ilevel
                   cell_igrid(ic) = igrid
                   cell_i1(ic) = i1
                   cell_i2(ic) = i2
                   cell_i3(ic) = i3
                end do
             end do
          end do
       end do
    end do

  end subroutine preset_cell_id

  logical function equal_grid_cell(a,b)
    implicit none
    type(grid_cell), intent(in) :: a,b
    equal_grid_cell = a%ic == b%ic
  end function equal_grid_cell

  type(grid_cell) function new_grid_cell_5d(i1, i2, i3, ilevel, igrid, geo) result(cell)
    implicit none
    integer,intent(in) :: i1, i2, i3, ilevel, igrid
    type(grid_geometry_desc),intent(in), target :: geo
    type(grid_desc),pointer :: grid
    grid => geo%levels(ilevel)%grids(igrid)
    cell%ic = (grid%start_id-1) + (i3-1)*grid%n1*grid%n2 + (i2-1)*grid%n1 + i1
    cell%ilevel = ilevel
    cell%igrid = igrid
    cell%i1 = i1
    cell%i2 = i2
    cell%i3 = i3
  end function new_grid_cell_5d

  type(grid_cell) function new_grid_cell_1d(ic, geo) result(cell)
    implicit none
    integer,intent(in) :: ic
    type(grid_geometry_desc),intent(in) :: geo
    cell%ic = ic
    cell%ilevel = cell_ilevel(ic)
    cell%igrid = cell_igrid(ic)
    cell%i1 = cell_i1(ic)
    cell%i2 = cell_i2(ic)
    cell%i3 = cell_i3(ic)
  end function new_grid_cell_1d

end module type_grid_cell
