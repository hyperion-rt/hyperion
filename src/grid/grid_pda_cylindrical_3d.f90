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
    select case(wall)
    case(1)
       geometrical_factor = 2. * geo%w1(cell%i1)/(geo%w1(cell%i1) + geo%w1(cell%i1+1))
    case(2)
       geometrical_factor = 2. * geo%w1(cell%i1+1)/(geo%w1(cell%i1) + geo%w1(cell%i1+1))
    case default
       geometrical_factor = 1._dp
    end select
  end function geometrical_factor

  subroutine check_allowed_pda(do_pda)
    implicit none
    logical,intent(inout) :: do_pda(:)
    integer :: i1, i2, i3, ic
    do i2=1,geo%n2
       do i3=1,geo%n3
          ic = cell_id(1, i2, i3, geo)
          do_pda(ic) = .false.
          ic = cell_id(geo%n1, i2, i3, geo)
          do_pda(ic) = .false.
       end do
    end do
    do i1=1,geo%n1
       do i3=1,geo%n3
          ic = cell_id(i1, 1, i3, geo)
          do_pda(ic) = .false.
          ic = cell_id(i1, geo%n2, i3, geo)
          do_pda(ic) = .false.
       end do
    end do
  end subroutine check_allowed_pda

end module grid_pda_geometry
