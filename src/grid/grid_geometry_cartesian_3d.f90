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
  public :: next_cell
  public :: place_in_cell
  public :: in_correct_cell
  public :: random_position_cell
  public :: find_wall
  public :: distance_to_closest_wall

  logical :: debug = .false.

  public :: escaped
  interface escaped
     module procedure escaped_photon
     module procedure escaped_cell
  end interface escaped

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

  type(grid_cell) function next_cell(cell, direction, intersection)
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
    next_cell = new_grid_cell(i1, i2, i3, geo)
  end function next_cell

  logical function in_correct_cell(p)
    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: frac
    icell_actual = find_cell(p)
    if(p%on_wall) then
       select case(p%on_wall_id)
       case(1)
          frac = (p%r%x - geo%w1(p%icell%i1)) / (geo%w1(p%icell%i1+1) - geo%w1(p%icell%i1))
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(2)
          frac = (p%r%x - geo%w1(p%icell%i1+1)) / (geo%w1(p%icell%i1+1) - geo%w1(p%icell%i1))
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(3)
          frac = (p%r%y - geo%w2(p%icell%i2)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(4)
          frac = (p%r%y - geo%w2(p%icell%i2+1)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(5)
          frac = (p%r%z - geo%w3(p%icell%i3)) / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i2 == p%icell%i2
       case(6)
          frac = (p%r%z - geo%w3(p%icell%i3+1)) / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i2 == p%icell%i2
       case default
          call warn("in_correct_cell","invalid on_wall_id")
          in_correct_cell = .false.
       end select
       in_correct_cell = abs(frac) < 1.e-3_dp .and. in_correct_cell
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

  subroutine find_wall(p,radial,tmin,id_min)

    implicit none

    type(photon), intent(inout) :: p
    ! Position and direction

    logical,intent(in) :: radial

    integer,intent(out) :: id_min
    ! ID of next wall

    real(dp),intent(out)  :: tmin
    ! tmin to nearest wall

    real(dp) :: tx, ty, tz

    logical :: pos_vx, pos_vy, pos_vz

    ! Can store this in photon, because it only changes at each interaction
    pos_vx = p%v%x > 0._dp ! whether photon is moving in the +ve x direction
    pos_vy = p%v%y > 0._dp ! whether photon is moving in the +ve y direction
    pos_vz = p%v%z > 0._dp ! whether photon is moving in the +ve z direction

    ! Store inv_v to go faster
    if(pos_vx) then
       tx = ( geo%w1(p%icell%i1+1) - p%r%x ) / p%v%x
    else if(p%v%x < 0._dp) then
       tx = ( geo%w1(p%icell%i1)   - p%r%x ) / p%v%x
    else
       tx = huge(1._dp)
    end if

    if(pos_vy) then
       ty = ( geo%w2(p%icell%i2+1) - p%r%y ) / p%v%y
    else if(p%v%y < 0._dp) then
       ty = ( geo%w2(p%icell%i2)   - p%r%y ) / p%v%y
    else
       ty = huge(1._dp)
    end if

    if(pos_vz) then
       tz = ( geo%w3(p%icell%i3+1) - p%r%z ) / p%v%z
    else if(p%v%z < 0._dp) then
       tz = ( geo%w3(p%icell%i3)   - p%r%z ) / p%v%z
    else
       tz = huge(1._dp)
    end if

    ! Following is potential slowdown, in fact, could just test tmin after
    if(min(tx,ty,tz) .lt. 0._dp) call error("find_wall","negative t")

    ! Find the closest of the three walls. The following effectively
    ! finds the minimum of three values. A lot of code for such a
    ! small thing, but this runs much much faster than using a built
    ! in min function or any kind of loop. Each iteraction comprises only
    ! three if statements and two pointer assignements.

    if(tx.lt.tz) then
       if(tx.lt.ty) then
          if(pos_vx) then
             id_min = 2
          else
             id_min = 1
          end if
          tmin = tx
       else
          if(pos_vy) then
             id_min = 4
          else
             id_min = 3
          end if
          tmin = ty
       end if
    else
       if(tz.lt.ty) then
          if(pos_vz) then
             id_min = 6
          else
             id_min = 5
          end if
          tmin = tz
       else
          if(pos_vy) then
             id_min = 4
          else
             id_min = 3
          end if
          tmin = ty
       end if
    end if

  end subroutine find_wall

end module grid_geometry_specific
