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
  public :: setup_grid_geometry

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

  type(grid_geometry_desc),public,target :: geo

contains

  logical pure function equal_nulp(x, y, n)
    implicit none
    real(8),intent(in) :: x, y
    integer,intent(in) :: n
    if(x == y) then
       equal_nulp = .true.
    else
       equal_nulp = abs(x - y) <= n * spacing(max(x, y))
    end if
  end function equal_nulp

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    select case(idir)
    case(1)
       cell_width = geo%dw(cell%i1)
    case(2)
       cell_width = geo%dz(cell%i2)
    case(3)
       cell_width = geo%w(cell%i1) * geo%dphi(cell%i3)
    end select
  end function cell_width

  real(dp) function cell_area(cell, iface)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: iface
    select case(iface)
    case(1)
       cell_area = geo%w1(cell%i1) * geo%dz(cell%i2) * geo%dphi(cell%i3)
    case(2)
       cell_area = geo%w1(cell%i1+1) * geo%dz(cell%i2) * geo%dphi(cell%i3)
    case(3,4)
       cell_area = 0.5 * geo%dw2(cell%i1) * geo%dphi(cell%i3)
    case(5,6)
       cell_area = geo%dw(cell%i1) * geo%dz(cell%i2)
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

    if(trim(geo%type).ne.'cyl_pol') call error("setup_grid_geometry","grid is not cylindrical polar")

    if(main_process()) write(*,'(" [setup_grid_geometry] Reading cylindrical polar grid")')

    call mp_table_read_column_auto(group, 'walls_1', 'w', geo%w1)
    call mp_table_read_column_auto(group, 'walls_2', 'z', geo%w2)
    call mp_table_read_column_auto(group, 'walls_3', 'p', geo%w3)

    if(any(geo%w1 < 0.)) then
       call error("setup_grid_geometry","w walls should be positive")
    end if

    if(any(geo%w3 < 0.).or.any(geo%w3 > twopi)) then
       call error("setup_grid_geometry","phi walls should be between 0 and 2*pi")
    end if

    geo%n1 = size(geo%w1) - 1
    geo%n2 = size(geo%w2) - 1
    geo%n3 = size(geo%w3) - 1
    geo%n_cells = geo%n1 * geo%n2 * geo%n3
    geo%n_masked = geo%n_cells

    allocate(geo%w(geo%n1))
    allocate(geo%dw(geo%n1))
    allocate(geo%dw2(geo%n1))
    allocate(geo%dz(geo%n2))
    allocate(geo%dphi(geo%n3))

    where(geo%w1(:geo%n1) == 0.)
       geo%w = geo%w1(2:) / 2._dp
    elsewhere
       geo%w = 10._dp**((log10(geo%w1(:geo%n1)) + log10(geo%w1(2:))) / 2._dp)
    end where

    geo%dw   = geo%w1(2:)    - geo%w1(:geo%n1)
    geo%dw2  = geo%w1(2:)**2 - geo%w1(:geo%n1)**2
    geo%dz   = geo%w2(2:)    - geo%w2(:geo%n2)
    geo%dphi = geo%w3(2:)    - geo%w3(:geo%n3)

    allocate(geo%volume(geo%n_cells))

    ! Compute cell volumes
    do ic=1,geo%n_cells
       cell = new_grid_cell(ic, geo)
       geo%volume(ic) = geo%dw2(cell%i1) * geo%dz(cell%i2) * geo%dphi(cell%i3) / 2._dp
    end do

    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')
    if(any(geo%dw==0._dp)) call error('setup_grid_geometry','all dw values should be greater than zero')
    if(any(geo%dz==0._dp)) call error('setup_grid_geometry','all dz values should be greater than zero')
    if(any(geo%dphi==0._dp)) call error('setup_grid_geometry','all dphi values should be greater than zero')

    ! Compute other useful quantities

    allocate(geo%wr2(geo%n1+1))
    geo%wr2 = geo%w1**2.

    allocate(geo%wtanp(geo%n3+1))
    geo%wtanp = tan(geo%w3)

    if(geo%n3 == 1) then
       geo%n_dim = 2
    else
       geo%n_dim = 3
    end if

    allocate(geo%ew1(geo%n1 + 1))
    allocate(geo%ew2(geo%n2 + 1))
    allocate(geo%ew3(geo%n3 + 1))

    geo%ew1 = 3 * spacing(geo%w1)
    geo%ew2 = 3 * spacing(geo%w2)
    geo%ew3 = 3 * spacing(1._dp) ! just use 1 since angles are of that order

  end subroutine setup_grid_geometry

  subroutine grid_geometry_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_geometry_debug

  type(grid_cell) function find_cell(p) result(icell)

    implicit none

    type(photon),intent(in) :: p
    real(dp) :: w_sq,phi
    integer :: i1, i2, i3

    if(debug) write(*,'(" [debug] find_cell")')

    w_sq = p%r%x*p%r%x+p%r%y*p%r%y

    ! If the photon is at the grid origin, then there is no way to determine
    ! the correct phi cell from the position. However, the direction vector
    ! indicates the direction the photon will move in, and since the grid cell
    ! walls do not curve as seen from the origin, we can use the direction
    ! vector to find out the correct cell.
    if(w_sq == 0._dp) then
       phi = atan2(p%v%y,p%v%x)
       if(phi < 0._dp) phi = phi + twopi
    else
       phi = atan2(p%r%y,p%r%x)
       if(phi < 0._dp) phi = phi + twopi
    end if

    ! Find the cells that the photon falls in
    i1 = locate(geo%wr2,w_sq)
    i2 = locate(geo%w2,p%r%z)
    i3 = locate(geo%w3,phi)

    ! We now check whether the photon falls outside the grid in any of the
    ! dimensions.

    if(i1<1.or.i1>geo%n1) then
       call warn("find_cell","photon not in grid (in r direction)")
       icell = invalid_cell
       return
    end if

    if(i2<1.or.i2>geo%n2) then
       call warn("find_cell","photon not in grid (in z direction)")
       icell = invalid_cell
       return
    end if

    if(i3<1.or.i3>geo%n3) then
       call warn("find_cell","photon not in grid (in phi direction)")
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

    real(dp) :: w_sq
    real(dp) :: phi, phi_v, dphi
    integer,parameter :: eps = 3

    ! Initialize values. We start off by assuming the photon is not on any
    ! walls.
    p%on_wall = .false.
    p%on_wall_id = no_wall

    ! Pre-compute useful quantities
    w_sq = p%r%x*p%r%x+p%r%y*p%r%y

    if(w_sq == 0._dp) then
       phi = atan2(p%v%y,p%v%x)
       if(phi < 0._dp) phi = phi + twopi
    else
       phi = atan2(p%r%y,p%r%x)
       if(phi < 0._dp) phi = phi + twopi
    end if

    ! Find whether the photon is on a radial wall.

    if((p%r%x * p%v%x + p%r%y * p%v%y) >= 0._dp) then
       if(equal_nulp(w_sq, geo%wr2(p%icell%i1), eps)) then ! photon is on inner wall
          p%on_wall_id%w1 = -1
       else if(equal_nulp(w_sq, geo%wr2(p%icell%i1 + 1), eps)) then ! photon is on outer wall
          p%on_wall_id%w1 = -1
          p%icell%i1 = p%icell%i1 + 1
       end if
    else
       if(equal_nulp(w_sq, geo%wr2(p%icell%i1), eps)) then ! photon is on inner wall
          p%on_wall_id%w1 = +1
          p%icell%i1 = p%icell%i1 - 1
       else if(equal_nulp(w_sq, geo%wr2(p%icell%i1 + 1), eps)) then ! photon is on outer wall
          p%on_wall_id%w1 = +1
       end if
    end if

    ! Find whether the photon is on a vertical wall

    if(p%v%z > 0._dp) then
       if(equal_nulp(p%r%z, geo%w2(p%icell%i2), eps)) then ! photon is on lower wall
          p%on_wall_id%w2 = -1
       else if(equal_nulp(p%r%z, geo%w2(p%icell%i2 + 1), eps)) then ! photon is on upper wall
          p%on_wall_id%w2 = -1
          p%icell%i2 = p%icell%i2 + 1
       end if
    else if(p%v%z < 0._dp) then
       if(equal_nulp(p%r%z, geo%w2(p%icell%i2), eps)) then ! photon is on lower wall
          p%on_wall_id%w2 = +1
          p%icell%i2 = p%icell%i2 - 1
       else if(equal_nulp(p%r%z, geo%w2(p%icell%i2 + 1), eps)) then ! photon is on upper wall
          p%on_wall_id%w2 = +1
       end if
    end if

    ! Find whether the photon is on an azimuthal wall

    if(p%r%x == 0._dp .and. p%r%y == 0._dp .and. p%v%x == 0._dp .and. p%v%y == 0._dp) then

       ! don't place on wall as on origin and moving up, so on *all* phi walls simultaneously

    else if(equal_nulp(phi, geo%w3(p%icell%i3), eps)) then  ! photon is on inner wall

       ! Find the angle between the direction of motion and the wall.
       phi_v = atan2(p%v%y, p%v%x)
       dphi = phi_v - geo%w3(p%icell%i3)
       if(dphi < -pi) dphi = dphi + twopi

       if(dphi > 0._dp) then  ! photon is moving towards upper cells
          p%on_wall_id%w3 = -1
       else  ! photon is moving towards lower cells
          p%on_wall_id%w3 = +1
          p%icell%i3 = p%icell%i3 - 1
          if(p%icell%i3==0) p%icell%i3 = geo%n3
       end if

    else if(equal_nulp(phi, geo%w3(p%icell%i3 + 1), eps)) then  ! photon is on outer wall

       ! Find the angle between the direction of motion and the wall.
       phi_v = atan2(p%v%y, p%v%x)
       dphi = phi_v - geo%w3(p%icell%i3 + 1)
       if(dphi < -pi) dphi = dphi + twopi

       if(dphi > 0._dp) then  ! photon is moving towards upper cells
          p%on_wall_id%w3 = -1
          p%icell%i3 = p%icell%i3 + 1
          if(p%icell%i3==geo%n3+1) p%icell%i3 = 1
       else  ! photon is moving towards lower cells
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
    if(cell%i1<1) return
    if(cell%i1>geo%n1) return
    if(cell%i2<1) return
    if(cell%i2>geo%n2) return
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
       if(i3==0) i3 = geo%n3
    case(6)
       i3 = i3 + 1
       if(i3==geo%n3+1) i3 = 1
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
       if(i3==0) i3 = geo%n3
    else if(direction%w3 == +1) then
       i3 = i3 + 1
       if(i3==geo%n3+1) i3 = 1
    end if
    next_cell_wall_id = new_grid_cell(i1, i2, i3, geo)
  end function next_cell_wall_id

  logical function in_correct_cell(p)

    implicit none

    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: phi,frac,dphi,r_sq,w_sq
    real(dp),parameter :: threshold = 1.e-3_dp

    icell_actual = find_cell(p)

    if(p%on_wall) then

       in_correct_cell = .true.

       w_sq = p%r%x*p%r%x+p%r%y*p%r%y

       if(w_sq == 0._dp) then
          phi = atan2(p%v%y,p%v%x)
          if(phi < 0._dp) phi = phi + twopi
       else
          phi = atan2(p%r%y,p%r%x)
          if(phi < 0._dp) phi = phi + twopi
       end if

       if(p%on_wall_id%w1 == -1) then
          if(geo%w1(p%icell%i1) .ne. sqrt(w_sq)) then
             frac = sqrt(w_sq) / geo%w1(p%icell%i1) - 1._dp
             in_correct_cell = in_correct_cell .and. abs(frac) < threshold
          end if
       else if(p%on_wall_id%w1 == +1) then
          if(geo%w1(p%icell%i1 + 1) .ne. sqrt(w_sq)) then
             frac = sqrt(w_sq) / geo%w1(p%icell%i1+1) - 1._dp
             in_correct_cell = in_correct_cell .and. abs(frac) < threshold
          end if
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i1 == p%icell%i1
       end if

       if(p%on_wall_id%w2 == -1) then
          frac = (p%r%z - geo%w2(p%icell%i2)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w2 == +1) then
          frac = (p%r%z - geo%w2(p%icell%i2+1)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i2 == p%icell%i2
       end if

       if(p%on_wall_id%w3 == -1) then
          dphi = phi - geo%w3(p%icell%i3)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
          frac = dphi / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w3 == +1) then
          dphi = phi - geo%w3(p%icell%i3+1)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
          frac = dphi / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
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
    real(dp) :: r,z,p

    call random(r)
    call random(z)
    call random(p)
    r = sqrt(r * (geo%w1(icell%i1+1)**2. - geo%w1(icell%i1)**2.) + geo%w1(icell%i1)**2.)
    z = z * (geo%w2(icell%i2+1) - geo%w2(icell%i2)) + geo%w2(icell%i2)
    p = p * (geo%w3(icell%i3+1) - geo%w3(icell%i3)) + geo%w3(icell%i3)

    ! If the photon doesn't fall within the walls, the cell is probably very
    ! narrow along that direction, so can just choose the mid-point

    if(r <= geo%w1(icell%i1) .or. r >= geo%w1(icell%i1+1)) then
       r = 0.5_dp * (geo%w1(icell%i1) + geo%w1(icell%i1+1))
    end if

    if(z <= geo%w2(icell%i2) .or. z >= geo%w2(icell%i2+1)) then
       z = 0.5_dp * (geo%w2(icell%i2) + geo%w2(icell%i2+1))
    end if

    if(p <= geo%w3(icell%i3) .or. p >= geo%w3(icell%i3+1)) then
       p = 0.5_dp * (geo%w3(icell%i3) + geo%w3(icell%i3+1))
    end if

    pos%x = r * cos(p)
    pos%y = r * sin(p)
    pos%z = z

  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p
    real(dp) :: r
    real(dp) :: d1,d2,d3,d4,d5,d6

    ! Cylindrical walls - point-cylinder distance

    r = sqrt(p%r%x*p%r%x+p%r%y*p%r%y)
    d1 = r - geo%w1(p%icell%i1)
    d2 = geo%w1(p%icell%i1+1) - r

    ! z walls - point-plane distance

    d3 = p%r%z - geo%w2(p%icell%i2)
    d4 = geo%w2(p%icell%i2+1) - p%r%z

    ! Phi walls - point-plane distance

    if(geo%n_dim == 3) then
       d5 = abs(geo%wtanp(p%icell%i3) * p%r%x - p%r%y) / sqrt(geo%wtanp(p%icell%i3)**2 + 1._dp)
       d6 = abs(geo%wtanp(p%icell%i3+1) * p%r%x - p%r%y) / sqrt(geo%wtanp(p%icell%i3+1)**2 + 1._dp)
    else
       d5 = huge(1._dp)
       d6 = huge(1._dp)
    end if

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

    real(dp) :: pB,pC,pC_1,pC_2,t1,t2
    real(dp) :: v2_xy
    real(dp) :: rv_xy
    real(dp) :: r2_xy

    real(dp) :: x_i,y_i,phi_i, dphi

    ! Find the intersections with all walls:

    call reset_t()

    ! -------------------------------------------------
    ! Compute general values
    ! -------------------------------------------------

    v2_xy = p%v%x*p%v%x+p%v%y*p%v%y

    rv_xy = p%r%x*p%v%x+p%r%y*p%v%y

    r2_xy = p%r%x*p%r%x+p%r%y*p%r%y

    ! -------------------------------------------------
    ! Cylindrical walls
    ! -------------------------------------------------

    ! Pre-compute some coefficients of the quadratic that don't depend on
    ! the wall, only the position and direction of the photon.
    pB=rv_xy / v2_xy ; pB = pB + pB
    pC=r2_xy / v2_xy

    ! Check for intersections with inner wall.
    pC_1 = pC - geo%wr2(p%icell%i1) / v2_xy
    call quadratic_pascal_reduced(pB,pC_1,t1,t2)

    ! If we are on the inner wall, then we should discard the
    ! intersection with the smallest absolute value as this will be the
    ! wall we are on. Otherwise, we should include both values.
    if(p%on_wall_id%w1 == -1) then
       if(abs(t1) < abs(t2)) then
          call insert_t(t2, 1, -1, geo%ew1(p%icell%i1))
       else
          call insert_t(t1, 1, -1, geo%ew1(p%icell%i1))
       end if
    else
       call insert_t(t1,1, -1,geo%ew1(p%icell%i1))
       call insert_t(t2,1, -1,geo%ew1(p%icell%i1))
    end if

    ! Check for intersections with outer wall.
    pC_2 = pC - geo%wr2(p%icell%i1+1) / v2_xy
    call quadratic_pascal_reduced(pB,pC_2,t1,t2)

    ! If we are on the outer wall, then we should discard the
    ! intersection with the smallest absolute value as this will be the
    ! wall we are on. Otherwise, we should include both values.
    if(p%on_wall_id%w1 == +1) then
       if(abs(t1) < abs(t2)) then
          call insert_t(t2, 1, +1, geo%ew1(p%icell%i1 + 1))
       else
          call insert_t(t1, 1, +1, geo%ew1(p%icell%i1 + 1))
       end if
    else
       call insert_t(t1,1, +1, geo%ew1(p%icell%i1 + 1))
       call insert_t(t2,1, +1, geo%ew1(p%icell%i1 + 1))
    end if

    ! -------------------------------------------------
    ! z walls
    ! -------------------------------------------------

    ! Check for intersection with lower wall
    if(p%on_wall_id%w2 /= -1) then
       t1 = ( geo%w2(p%icell%i2)   - p%r%z ) / p%v%z
       call insert_t(t1,2, -1, 0._dp)
    end if

    ! Check for intersection with upper wall
    if(p%on_wall_id%w2 /= +1) then
       t2 = ( geo%w2(p%icell%i2+1) - p%r%z ) / p%v%z
       call insert_t(t2,2, +1, 0._dp)
    end if

    ! -------------------------------------------------
    ! phi walls
    ! -------------------------------------------------

    ! For performance reasons, we only check for intersections with phi
    ! walls if there are multiple cells in phi
    if(geo%n_dim == 3) then

       ! If the photon is on a phi wall and moving in the direction of the
       ! wall, then we don't want to check for intersections, but we want to
       ! make sure the photon is still considered to be on the wall
       ! afterwards. Therefore, we calculate dphi to find the angle between
       ! the direction of motion and the wall the photon is on.

       if(p%on_wall_id%w3 == -1) then
          dphi = atan2(p%v%y, p%v%x) - geo%w3(p%icell%i3)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
       end if
       if(p%on_wall_id%w3 == +1) then
          dphi = atan2(p%v%y, p%v%x) - geo%w3(p%icell%i3+1)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
       end if

       if(p%on_wall_id%w3 == +1 .and. abs(dphi) < geo%ew3(p%icell%i3+1)) then

          iext%w3 = +1

       else if(p%on_wall_id%w3 == -1 .and. abs(dphi) < geo%ew3(p%icell%i3)) then

          iext%w3 = -1

       ! Only check for intersections if the current position is not along
       ! the (x, y) = (0, 0) axis. If the photon is along this axis, then
       ! either it is and will remain on a wall until the next iteraction, or
       ! it is not and will not intersect a phi wall.
       else if(r2_xy > 0._dp) then

          ! Only check for intersections if we are not on the wall. We can
          ! do this for the phi walls because they are planes, so if we are on
          ! them, we can't intersect them.
          if(p%on_wall_id%w3 /= -1) then

             ! Find intersection with lower phi wall
             t1 = - ( geo%wtanp(p%icell%i3) * p%r%x - p%r%y ) / ( geo%wtanp(p%icell%i3) * p%v%x - p%v%y )

             ! Find position of intersection in x,y
             x_i = p%r%x + p%v%x * t1
             y_i = p%r%y + p%v%y * t1
             phi_i = atan2(y_i, x_i)
             dphi = abs(phi_i - geo%w3(p%icell%i3))
             if(dphi > pi) dphi = abs(dphi - twopi)
             if(dphi < 0.5 * pi) call insert_t(t1, 3, -1, 0._dp)

          end if

          ! Only check for intersections if we are not on the wall. We can
          ! do this for the phi walls because they are planes, so if we are on
          ! them, we can't intersect them.
          if(p%on_wall_id%w3 /= +1) then

             ! Find intersection with upper phi wall
             t2 = - ( geo%wtanp(p%icell%i3+1) * p%r%x - p%r%y ) / ( geo%wtanp(p%icell%i3+1) * p%v%x - p%v%y )

             ! Find position of intersection in x,y
             x_i = p%r%x + p%v%x * t2
             y_i = p%r%y + p%v%y * t2
             phi_i = atan2(y_i, x_i)
             dphi = abs(phi_i - geo%w3(p%icell%i3+1))
             if(dphi > pi) dphi = abs(dphi - twopi)
             if(dphi < 0.5 * pi) call insert_t(t2, 3, +1, 0._dp)

          end if

       end if

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
