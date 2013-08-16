module grid_geometry_specific

  use core_lib
  use mpi_core
  use mpi_hdf5_io
  use type_photon
  use type_grid_cell
  use type_grid
  use grid_io
  use counters

  implicit none
  save

  private

  ! Photon position routines
  public :: cell_width
  public :: grid_geometry_debug
  public :: find_cell
  public :: place_in_cell
  public :: in_correct_cell
  public :: random_position_cell
  public :: find_wall
  public :: distance_to_closest_wall

  logical :: debug = .false.

  real(dp) :: tmin, emin
  type(wall_id) :: imin, iext

  public :: escaped
  interface escaped
     module procedure escaped_photon
     module procedure escaped_cell
  end interface escaped

  public :: next_cell
  interface next_cell
     module procedure next_cell_int
     module procedure next_cell_wall_id
  end interface next_cell

  public :: setup_grid_geometry
  type(grid_geometry_desc),public,target :: geo

contains

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    select case(idir)
    case(1)
       cell_width = geo%dx(cell%i1)
    case(2)
       cell_width = geo%dy(cell%i2)
    case(3)
       cell_width = geo%dz(cell%i3)
    end select
  end function cell_width

  real(dp) function cell_area(cell, iface)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: iface
    select case(iface)
    case(1,2)
       cell_area = geo%dy(cell%i2) * geo%dz(cell%i3)
    case(3,4)
       cell_area = geo%dz(cell%i3) * geo%dx(cell%i1)
    case(5,6)
       cell_area = geo%dx(cell%i1) * geo%dy(cell%i2)
    end select
  end function cell_area

  subroutine setup_grid_geometry(group)

    implicit none

    integer(hid_t),intent(in) :: group
    integer :: ic
    type(grid_cell) :: cell

    ! Read geometry file
    call mp_read_keyword(group, '.', "geometry", geo%id)
    call mp_read_keyword(group, '.', "grid_type", geo%type)

    if(trim(geo%type).ne.'car') call error("setup_grid_geometry","grid is not cartesian")

    if(main_process()) write(*,'(" [setup_grid_geometry] Reading cartesian grid")')

    call mp_table_read_column_auto(group, 'walls_1', 'x', geo%w1)
    call mp_table_read_column_auto(group, 'walls_2', 'y', geo%w2)
    call mp_table_read_column_auto(group, 'walls_3', 'z', geo%w3)

    geo%n1 = size(geo%w1) - 1
    geo%n2 = size(geo%w2) - 1
    geo%n3 = size(geo%w3) - 1
    geo%n_cells = geo%n1 * geo%n2 * geo%n3
    geo%n_masked = geo%n_cells

    allocate(geo%dx(geo%n1))
    allocate(geo%dy(geo%n2))
    allocate(geo%dz(geo%n3))

    geo%dx = geo%w1(2:) - geo%w1(:geo%n1)
    geo%dy = geo%w2(2:) - geo%w2(:geo%n2)
    geo%dz = geo%w3(2:) - geo%w3(:geo%n3)

    allocate(geo%volume(geo%n_cells))

    ! Compute cell volumes
    do ic=1,geo%n_cells
       cell = new_grid_cell(ic, geo)
       geo%volume(ic) = geo%dx(cell%i1) * geo%dy(cell%i2) * geo%dz(cell%i3)
    end do

    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')
    if(any(geo%dx==0._dp)) call error('setup_grid_geometry','all dx values should be greater than zero')
    if(any(geo%dy==0._dp)) call error('setup_grid_geometry','all dy values should be greater than zero')
    if(any(geo%dz==0._dp)) call error('setup_grid_geometry','all dz values should be greater than zero')

    ! Compute other useful quantities
    geo%n_dim = 3

    allocate(geo%ew1(geo%n1 + 1))
    allocate(geo%ew2(geo%n2 + 1))
    allocate(geo%ew3(geo%n3 + 1))

    geo%ew1 = 3 * spacing(geo%w1)
    geo%ew2 = 3 * spacing(geo%w2)
    geo%ew3 = 3 * spacing(geo%w3)

  end subroutine setup_grid_geometry

  subroutine grid_geometry_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_geometry_debug

  type(grid_cell) function find_cell(p) result(icell)
    implicit none
    type(photon),intent(in) :: p
    integer :: i1, i2, i3
    if(debug) write(*,'(" [debug] find_cell")')
    i1 = locate(geo%w1,p%r%x)
    i2 = locate(geo%w2,p%r%y)
    i3 = locate(geo%w3,p%r%z)
    if(i1<1.or.i1>geo%n1) then
       call warn("find_cell","photon not in cell (in x direction)")
       icell = invalid_cell
       return
    end if
    if(i2<1.or.i2>geo%n2) then
       call warn("find_cell","photon not in cell (in y direction)")
       icell = invalid_cell
       return
    end if
    if(i3<1.or.i3>geo%n3) then
       call warn("find_cell","photon not in cell (in z direction)")
       icell = invalid_cell
       return
    end if
    icell = new_grid_cell(i1, i2, i3, geo)
  end function find_cell

  subroutine adjust_wall(p)

    ! In future, this could be called at peeloff time instead of
    ! place_inside_cell, but if we want to do that, we need this subroutine to
    ! use locate to find the correct cell if the velocity is zero along one of
    ! the components, so that it is reset to the 'default' find_cell value.

    implicit none

    type(photon), intent(inout) :: p

    ! Initialize values
    p%on_wall = .false.
    p%on_wall_id = no_wall

    ! Find whether the photon is on an x-wall
    if(p%v%x > 0._dp) then
       if(p%r%x == geo%w1(p%icell%i1)) then
          p%on_wall_id%w1 = -1
       else if(p%r%x == geo%w1(p%icell%i1 + 1)) then
          p%on_wall_id%w1 = -1
          p%icell%i1 = p%icell%i1 + 1
       end if
    else if(p%v%x < 0._dp) then
       if(p%r%x == geo%w1(p%icell%i1)) then
          p%on_wall_id%w1 = +1
          p%icell%i1 = p%icell%i1 - 1
       else if(p%r%x == geo%w1(p%icell%i1 + 1)) then
          p%on_wall_id%w1 = +1
       end if
    end if

    ! Find whether the photon is on a y-wall
    if(p%v%y > 0._dp) then
       if(p%r%y == geo%w2(p%icell%i2)) then
          p%on_wall_id%w2 = -1
       else if(p%r%y == geo%w2(p%icell%i2 + 1)) then
          p%on_wall_id%w2 = -1
          p%icell%i2 = p%icell%i2 + 1
       end if
    else if(p%v%y < 0._dp) then
       if(p%r%y == geo%w2(p%icell%i2)) then
          p%on_wall_id%w2 = +1
          p%icell%i2 = p%icell%i2 - 1
       else if(p%r%y == geo%w2(p%icell%i2 + 1)) then
          p%on_wall_id%w2 = +1
       end if
    end if

    ! Find whether the photon is on a z-wall
    if(p%v%z > 0._dp) then
       if(p%r%z == geo%w3(p%icell%i3)) then
          p%on_wall_id%w3 = -1
       else if(p%r%z == geo%w3(p%icell%i3 + 1)) then
          p%on_wall_id%w3 = -1
          p%icell%i3 = p%icell%i3 + 1
       end if
    else if(p%v%z < 0._dp) then
       if(p%r%z == geo%w3(p%icell%i3)) then
          p%on_wall_id%w3 = +1
          p%icell%i3 = p%icell%i3 - 1
       else if(p%r%z == geo%w3(p%icell%i3 + 1)) then
          p%on_wall_id%w3 = +1
       end if
    end if

    p%on_wall = p%on_wall_id%w1 /= 0 .or. p%on_wall_id%w2 /= 0 .or. p%on_wall_id%w3 /= 0

  end subroutine adjust_wall

  subroutine place_in_cell(p)

    implicit none

    type(photon),intent(inout) :: p

    p%icell = find_cell(p)

    if(p%icell == invalid_cell) then
       call warn("place_in_cell","place_in_cell failed - killing")
       killed_photons_geo = killed_photons_geo + 1
       p%killed = .true.
    else
       p%in_cell = .true.
    end if

    call adjust_wall(p)

  end subroutine place_in_cell

  logical function escaped_photon(p)
    implicit none
    type(photon),intent(in) :: p
    escaped_photon = escaped_cell(p%icell)
  end function escaped_photon

  logical function escaped_cell(cell)
    implicit none
    type(grid_cell),intent(in) :: cell
    escaped_cell = .true.
    if(cell%i1 < 1 .or. cell%i1 > geo%n1) return
    if(cell%i2 < 1 .or. cell%i2 > geo%n2) return
    if(cell%i3 < 1 .or. cell%i3 > geo%n3) return
    escaped_cell = .false.
  end function escaped_cell

  type(grid_cell) function next_cell_int(cell, direction, intersection)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    integer :: i1, i2, i3
    i1 = cell%i1
    i2 = cell%i2
    i3 = cell%i3
    select case(direction)
    case(1)
       i1 = i1 - 1
    case(2)
       i1 = i1 + 1
    case(3)
       i2 = i2 - 1
    case(4)
       i2 = i2 + 1
    case(5)
       i3 = i3 - 1
    case(6)
       i3 = i3 + 1
    end select
    next_cell_int = new_grid_cell(i1, i2, i3, geo)
  end function next_cell_int

  type(grid_cell) function next_cell_wall_id(cell, direction, intersection)
    implicit none
    type(grid_cell),intent(in) :: cell
    type(wall_id),intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    integer :: i1, i2, i3
    i1 = cell%i1
    i2 = cell%i2
    i3 = cell%i3
    if(direction%w1 == -1) then
       i1 = i1 - 1
    else if(direction%w1 == +1) then
       i1 = i1 + 1
    end if
    if(direction%w2 == -1) then
       i2 = i2 - 1
    else if(direction%w2 == +1) then
       i2 = i2 + 1
    end if
    if(direction%w3 == -1) then
       i3 = i3 - 1
    else if(direction%w3 == +1) then
       i3 = i3 + 1
    end if
    next_cell_wall_id = new_grid_cell(i1, i2, i3, geo)
  end function next_cell_wall_id

  logical function in_correct_cell(p)

    implicit none

    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: frac
    real(dp),parameter :: threshold = 1.e-3_dp

    icell_actual = find_cell(p)

    if(p%on_wall) then

       in_correct_cell = .true.

       if(p%on_wall_id%w1 == -1) then
          frac = (p%r%x - geo%w1(p%icell%i1)) / (geo%w1(p%icell%i1+1) - geo%w1(p%icell%i1))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w1 == +1) then
          frac = (p%r%x - geo%w1(p%icell%i1+1)) / (geo%w1(p%icell%i1+1) - geo%w1(p%icell%i1))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i1 == p%icell%i1
       end if

       if(p%on_wall_id%w2 == -1) then
          frac = (p%r%y - geo%w2(p%icell%i2)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w2 == +1) then
          frac = (p%r%y - geo%w2(p%icell%i2+1)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i2 == p%icell%i2
       end if

       if(p%on_wall_id%w3 == -1) then
          frac = (p%r%z - geo%w3(p%icell%i3)) / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w3 == +1) then
          frac = (p%r%z - geo%w3(p%icell%i3+1)) / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i3 == p%icell%i3
       end if

    else

       in_correct_cell = icell_actual == p%icell

    end if

  end function in_correct_cell

  subroutine random_position_cell(icell,pos)
    implicit none
    type(grid_cell),intent(in) :: icell
    type(vector3d_dp), intent(out) :: pos
    real(dp) :: x,y,z
    call random(x)
    call random(y)
    call random(z)
    pos%x = x * (geo%w1(icell%i1+1) - geo%w1(icell%i1)) + geo%w1(icell%i1)
    pos%y = y * (geo%w2(icell%i2+1) - geo%w2(icell%i2)) + geo%w2(icell%i2)
    pos%z = z * (geo%w3(icell%i3+1) - geo%w3(icell%i3)) + geo%w3(icell%i3)
  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p

    real(dp) :: d1,d2,d3,d4,d5,d6

    d1 = p%r%x - geo%w1(p%icell%i1)
    d2 = geo%w1(p%icell%i1+1) - p%r%x
    d3 = p%r%y - geo%w2(p%icell%i2)
    d4 = geo%w2(p%icell%i2+1) - p%r%y
    d5 = p%r%z - geo%w3(p%icell%i3)
    d6 = geo%w3(p%icell%i3+1) - p%r%z

    ! Find the smallest of the distances

    d = min(d1,d2,d3,d4,d5,d6)

    ! The closest distance should never be negative

    if(d < 0._dp) then
       call warn("distance_to_closest_wall","distance to closest wall is negative (assuming zero)")
       d = 0._dp
    end if

  end function distance_to_closest_wall

  subroutine find_wall(p,radial,tnearest,id_min)

    implicit none

    type(photon), intent(inout) :: p
    ! Position and direction

    logical,intent(in) :: radial

    type(wall_id),intent(out) :: id_min
    ! ID of next wall

    real(dp),intent(out)  :: tnearest
    ! tmin to nearest wall

    real(dp) :: t1, t2

    logical :: pos_vx, pos_vy, pos_vz

    call reset_t()

    if(p%on_wall_id%w1 /= -1) then
       t1 = ( geo%w1(p%icell%i1) - p%r%x ) / p%v%x
       call insert_t(t1,1, -1, geo%ew1(p%icell%i1))
    end if
    if(p%on_wall_id%w1 /= +1) then
       t2 = ( geo%w1(p%icell%i1+1) - p%r%x ) / p%v%x
       call insert_t(t2,1, +1, geo%ew1(p%icell%i1 + 1))
    end if

    if(p%on_wall_id%w2 /= -1) then
       t1 = ( geo%w2(p%icell%i2) - p%r%y ) / p%v%y
       call insert_t(t1,2, -1, geo%ew2(p%icell%i2))
    end if
    if(p%on_wall_id%w2 /= +1) then
       t2 = ( geo%w2(p%icell%i2+1) - p%r%y ) / p%v%y
       call insert_t(t2,2, +1, geo%ew2(p%icell%i2 + 1))
    end if

    if(p%on_wall_id%w3 /= -1) then
       t1 = ( geo%w3(p%icell%i3) - p%r%z ) / p%v%z
       call insert_t(t1,3, -1, geo%ew3(p%icell%i3))
    end if
    if(p%on_wall_id%w3 /= +1) then
       t2 = ( geo%w3(p%icell%i3+1) - p%r%z ) / p%v%z
       call insert_t(t2,3, +1, geo%ew3(p%icell%i3 + 1))
    end if

    call find_next_wall(tnearest,id_min)

  end subroutine find_wall

  subroutine reset_t()
    implicit none
    tmin = +huge(tmin)
    emin = 0.
    imin = no_wall
    iext = no_wall
  end subroutine reset_t

  subroutine insert_t(t, iw, i, e)
    implicit none
    real(dp),intent(in)    :: t, e
    integer,intent(in)    :: i, iw
    real(dp) :: emax
    if(debug) print *,'[debug] inserting t,i=',t, e, iw, i
    if(t > 0._dp) then
       emax = max(e, emin)
       if(t < tmin - emax) then
          tmin = t
          imin = no_wall
          emin = emax
          if(iw == 1) then
             imin%w1 = i
          else if(iw == 2) then
             imin%w2 = i
          else
             imin%w3 = i
          end if
       else if(t < tmin + emax) then
          emin = emax
          if(iw == 1) then
             imin%w1 = i
          else if(iw == 2) then
             imin%w2 = i
          else
             imin%w3 = i
          end if
       end if
    end if
  end subroutine insert_t

  subroutine find_next_wall(t,i)
    implicit none
    real(dp),intent(out)    :: t
    type(wall_id),intent(out)    :: i
    t = tmin
    i = imin + iext
    if(debug) print *,'[debug] selecting t,i=',t,i
  end subroutine find_next_wall

end module grid_geometry_specific
