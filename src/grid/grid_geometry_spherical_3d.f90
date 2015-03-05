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
       cell_width = geo%dr(cell%i1)
    case(2)
       cell_width = geo%r(cell%i1) * geo%dt(cell%i2)
    case(3)
       cell_width = geo%r(cell%i1) * sin(geo%t(cell%i2)) * geo%dphi(cell%i3)
    end select
  end function cell_width

  real(dp) function cell_area(cell, iface)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: iface
    select case(iface)
    case(1)
       cell_area = geo%w1(cell%i1)**2 * geo%dcost(cell%i2) * geo%dphi(cell%i3)
    case(2)
       cell_area = geo%w1(cell%i1+1)**2 * geo%dcost(cell%i2) * geo%dphi(cell%i3)
    case(3)
       cell_area = 0.5_dp * geo%dr2(cell%i1) * sin(geo%w2(cell%i2)) * geo%dphi(cell%i3)
    case(4)
       cell_area = 0.5_dp * geo%dr2(cell%i1) * sin(geo%w2(cell%i2 + 1)) * geo%dphi(cell%i3)
    case(5,6)
       cell_area = 0.5_dp * geo%dr2(cell%i1) * geo%dt(cell%i2)
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

    if(trim(geo%type).ne.'sph_pol') call error("setup_grid_geometry","grid is not spherical polar")

    if(main_process()) write(*,'(" [setup_grid_geometry] Reading spherical polar grid")')

    call mp_table_read_column_auto(group, 'walls_1', 'r', geo%w1)
    call mp_table_read_column_auto(group, 'walls_2', 't', geo%w2)
    call mp_table_read_column_auto(group, 'walls_3', 'p', geo%w3)

    if(any(geo%w1 < 0.)) then
       call error("setup_grid_geometry","r walls should be positive")
    end if

    if(any(geo%w2 < 0.).or.any(geo%w2 > pi)) then
       call error("setup_grid_geometry","theta walls should be between 0 and pi")
    end if

    if(any(geo%w3 < 0.).or.any(geo%w3 > twopi)) then
       call error("setup_grid_geometry","phi walls should be between 0 and 2*pi")
    end if

    geo%n1 = size(geo%w1) - 1
    geo%n2 = size(geo%w2) - 1
    geo%n3 = size(geo%w3) - 1
    geo%n_cells = geo%n1 * geo%n2 * geo%n3
    geo%n_masked = geo%n_cells

    allocate(geo%r(geo%n1))
    allocate(geo%dr(geo%n1))
    allocate(geo%dr2(geo%n1))
    allocate(geo%dr3(geo%n1))
    allocate(geo%t(geo%n2))
    allocate(geo%dt(geo%n2))
    allocate(geo%dcost(geo%n2))
    allocate(geo%dphi(geo%n3))

    where(geo%w1(:geo%n1) == 0.)
       geo%r = geo%w1(2:) / 2._dp
    elsewhere
       geo%r = 10._dp**((log10(geo%w1(:geo%n1)) + log10(geo%w1(2:))) / 2._dp)
    end where

    geo%t = (geo%w2(:geo%n2) + geo%w2(2:)) / 2._dp

    geo%dr    = geo%w1(2:)           - geo%w1(:geo%n1)
    geo%dr2   = geo%w1(2:)**2        - geo%w1(:geo%n1)**2
    geo%dr3   = geo%w1(2:)**3        - geo%w1(:geo%n1)**3
    geo%dt    = geo%w2(2:)           - geo%w2(:geo%n2)
    geo%dcost = cos(geo%w2(:geo%n2)) - cos(geo%w2(2:))
    geo%dphi  = geo%w3(2:)           - geo%w3(:geo%n3)

    allocate(geo%volume(geo%n_cells))

    ! Compute cell volumes
    do ic=1,geo%n_cells
       cell = new_grid_cell(ic, geo)
       geo%volume(ic) = geo%dr3(cell%i1) * geo%dcost(cell%i2) * geo%dphi(cell%i3) / 3._dp
    end do

    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')
    if(any(geo%dr==0._dp)) call error('setup_grid_geometry','all dr values should be greater than zero')
    if(any(geo%dt==0._dp)) call error('setup_grid_geometry','all dt values should be greater than zero')
    if(any(geo%dphi==0._dp)) call error('setup_grid_geometry','all dphi values should be greater than zero')

    ! Compute other useful quantities

    allocate(geo%wr2(geo%n1+1))
    geo%wr2 = geo%w1*geo%w1

    allocate(geo%wtanp(geo%n3+1))
    geo%wtanp = tan(geo%w3)

    allocate(geo%wtant(geo%n2+1))
    geo%wtant = tan(geo%w2)

    allocate(geo%wcost(geo%n2+1))
    geo%wcost = cos(geo%w2)

    allocate(geo%wsint(geo%n2+1))
    geo%wsint = sin(geo%w2)

    allocate(geo%wtant2(geo%n2+1))
    geo%wtant2 = geo%wtant*geo%wtant

    ! Search for midplane wall
    if(any(abs(geo%w2-pi/2._dp)<1.e-6)) geo%midplane = minloc(abs(geo%w2-pi/2._dp),1)

    if(geo%n3 == 1) then
       geo%n_dim = 2
    else
       geo%n_dim = 3
    end if

    allocate(geo%ew1(geo%n1 + 1))
    allocate(geo%ew2(geo%n2 + 1))
    allocate(geo%ew3(geo%n3 + 1))

    geo%ew1 = 3 * spacing(geo%w1)
    geo%ew2 = 3 * spacing(1._dp) ! just use 1 since angles are of that order
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
    real(dp) :: w_sq, r_sq, theta, phi
    integer :: i1, i2, i3

    if(debug) write(*,'(" [debug] find_cell")')

    ! Pre-compute useful quantities
    r_sq = p%r.dot.p%r
    w_sq = p%r%x*p%r%x+p%r%y*p%r%y

    ! If the photon is at the grid origin, then there is no way to determine
    ! the correct theta cell from the position. However, the direction vector
    ! indicates the direction the photon will move in, and since the grid cell
    ! walls do not curve as seen from the origin, we can use the direction
    ! vector to find out the correct cell.
    if(r_sq == 0._dp) then
       theta = atan2(sqrt(p%v%x*p%v%x+p%v%y*p%v%y),p%v%z)
    else
       theta = atan2(sqrt(p%r%x*p%r%x+p%r%y*p%r%y),p%r%z)
    end if

    ! Similarly to theta, if the photon is at (x, y) = (0, 0) - regardless
    ! of the value of z - then the correct phi cell cannot be determined from
    ! the position, so we have to use the direction vector.
    if(w_sq == 0._dp) then
       phi = atan2(p%v%y,p%v%x)
       if(phi < 0._dp) phi = phi + twopi
    else
       phi = atan2(p%r%y,p%r%x)
       if(phi < 0._dp) phi = phi + twopi
    end if

    ! Find the cells that the photon falls in
    i1 = locate(geo%wr2,r_sq)
    i2 = locate(geo%w2,theta)
    i3 = locate(geo%w3,phi)

    ! We now check whether the photon falls outside the grid in any of the
    ! dimensions.

    if(i1<1.or.i1>geo%n1) then
       call warn("find_cell","photon not in grid (in r direction)")
       icell = invalid_cell
       return
    end if

    if(i2<1.or.i2>geo%n2) then
       call warn("find_cell","photon not in grid (in theta direction)")
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

    real(dp) :: r_sq, w_sq
    real(dp) :: theta, theta_v
    real(dp) :: phi, phi_v, dphi
    integer,parameter :: eps = 3

    ! Initialize values. We start off by assuming the photon is not on any
    ! walls.
    p%on_wall = .false.
    p%on_wall_id = no_wall

    ! Pre-compute useful quantities
    r_sq = p%r.dot.p%r
    w_sq = p%r%x*p%r%x+p%r%y*p%r%y

    if(r_sq == 0._dp) then
       theta = atan2(sqrt(p%v%x*p%v%x+p%v%y*p%v%y),p%v%z)
    else
       theta = atan2(sqrt(p%r%x*p%r%x+p%r%y*p%r%y),p%r%z)
    end if

    if(w_sq == 0._dp) then
       phi = atan2(p%v%y,p%v%x)
       if(phi < 0._dp) phi = phi + twopi
    else
       phi = atan2(p%r%y,p%r%x)
       if(phi < 0._dp) phi = phi + twopi
    end if

    ! Find whether the photon is on a radial wall.

    if((p%r .dot. p%v) >= 0._dp) then  ! photon is moving outwards
       if(equal_nulp(r_sq, geo%wr2(p%icell%i1), eps)) then ! photon is on inner wall
          p%on_wall_id%w1 = -1
       else if(equal_nulp(r_sq, geo%wr2(p%icell%i1 + 1), eps)) then ! photon is on outer wall
          p%on_wall_id%w1 = -1
          p%icell%i1 = p%icell%i1 + 1
       end if
    else  ! photon is moving inwards
       if(equal_nulp(r_sq, geo%wr2(p%icell%i1), eps)) then ! photon is on inner wall
          p%on_wall_id%w1 = +1
          p%icell%i1 = p%icell%i1 - 1
       else if(equal_nulp(r_sq, geo%wr2(p%icell%i1 + 1), eps)) then ! photon is on outer wall
          p%on_wall_id%w1 = +1
       end if
    end if

    ! Find whether the photon is on a theta wall.

    ! We first have to check whether we are at the origin, because in that
    ! case we cannot rely on the position of the photon to check which wall it
    ! is on.
    if(r_sq == 0._dp) then  ! photon is at origin

       ! We need to ensure that the photon is not moving straight up or
       ! straight down. If it is, then we do not do anything since we do not
       ! treat the theta=0 and theta=pi walls as real walls. If the magnitude
       ! of the vertical velocity is not 1, then the photon is not moving
       ! vertically (since v is normalized).
       if(abs(p%v%z) < 1._dp) then

          ! Find direction of motion of the photon.
          theta_v = atan2(sqrt(p%v%x*p%v%x + p%v%y*p%v%y), p%v%z)

          ! Compare to the position of the lower and upper walls, and set
          ! the wall ID accordingly. We don't need to worry about changing the
          ! cell of the photon, since it is moving along the wall, so it
          ! doesn't matter whether it is set to one cell or the other.
          if(equal_nulp(theta_v, geo%w2(p%icell%i2), eps)) then
             p%on_wall_id%w2 = -1
          else if(equal_nulp(theta_v, geo%w2(p%icell%i2 + 1), eps)) then
             p%on_wall_id%w2 = +1
          end if

       end if

    ! If the photon is not at the origin, then we check whether it is on the
    ! lower wall, but we only do that if the lower wall is not theta=0.
    else if(p%icell%i2 > 1 .and. equal_nulp(theta, geo%w2(p%icell%i2), eps)) then

       if(p%icell%i2 == geo%midplane) then  ! lower wall is grid midplane

          ! The lower wall is the midplane of the grid, so we can just use
          ! v_x to determine the direction of motion relative to the wall.
          if(p%v%z > 0._dp) then
             p%on_wall_id%w2 = +1
             p%icell%i2 = p%icell%i2 - 1
          else
             p%on_wall_id%w2 = -1
          end if

       else  ! lower wall is not grid midplane

          ! We need to do some vector arithmetic and dot products to
          ! determine whether the photon is moving towards the lower or upper
          ! cells.
          if((sqrt(w_sq) * p%v%z * geo%wtant(p%icell%i2) -  (p%r%x * p%v%x + p%r%y * p%v%y) < 0._dp) .eqv. p%r%z > 0._dp) then
             p%on_wall_id%w2 = -1
          else
             p%on_wall_id%w2 = +1
             p%icell%i2 = p%icell%i2 - 1
          end if

       end if

    ! If the photon is not at the origin, or on the lower wall, then we
    ! check whether it is on the upper wall, but we only do that if the lower
    ! wall is not theta=pi.
    else if(p%icell%i2 + 1 < geo%n2 + 1 .and. equal_nulp(theta, geo%w2(p%icell%i2 + 1), eps)) then

       if(p%icell%i2 + 1 == geo%midplane) then  ! upper wall is grid midplane

          ! The lower wall is the midplane of the grid, so we can just use
          ! v_x to determine the direction of motion relative to the wall.
          if(p%v%z > 0._dp) then
             p%on_wall_id%w2 = +1
          else
             p%on_wall_id%w2 = -1
             p%icell%i2 = p%icell%i2 + 1
          end if

       else  ! upper wall is not grid midplane

          ! We need to do some vector arithmetic and dot products to
          ! determine whether the photon is moving towards the lower or upper
          ! cells.
          if((sqrt(w_sq) * p%v%z * geo%wtant(p%icell%i2 + 1)  -  (p%r%x * p%v%x + p%r%y * p%v%y) < 0._dp) .eqv. p%r%z > 0._dp) then
             p%on_wall_id%w2 = -1
             p%icell%i2 = p%icell%i2 + 1
          else
             p%on_wall_id%w2 = +1
          end if
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
    if(cell%i1<1)   return
    if(cell%i1>geo%n1) return
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
    real(dp) :: theta,phi,frac,dphi,r_sq,w_sq
    real(dp),parameter :: threshold = 1.e-3_dp

    icell_actual = find_cell(p)

    if(p%on_wall) then

       in_correct_cell = .true.

       r_sq = p%r.dot.p%r

       ! If we are at the origin, then there isn't much point in checking,
       ! since this is difficult numerically, and the photon has not
       ! propagated anyway.
       if(r_sq == 0._dp) return

       w_sq = p%r%x*p%r%x+p%r%y*p%r%y

       if(r_sq == 0._dp) then
          theta = atan2(sqrt(p%v%x*p%v%x+p%v%y*p%v%y),p%v%z)
       else
          theta = atan2(sqrt(p%r%x*p%r%x+p%r%y*p%r%y),p%r%z)
       end if

       if(w_sq == 0._dp) then
          phi = atan2(p%v%y,p%v%x)
          if(phi < 0._dp) phi = phi + twopi
       else
          phi = atan2(p%r%y,p%r%x)
          if(phi < 0._dp) phi = phi + twopi
       end if

       if(p%on_wall_id%w1 == -1) then
          if(geo%w1(p%icell%i1) .ne. sqrt(r_sq)) then
             frac = sqrt(r_sq) / geo%w1(p%icell%i1) - 1._dp
             in_correct_cell = in_correct_cell .and. abs(frac) < threshold
          end if
       else if(p%on_wall_id%w1 == +1) then
          if(geo%w1(p%icell%i1 + 1) .ne. sqrt(r_sq)) then
             frac = sqrt(r_sq) / geo%w1(p%icell%i1+1) - 1._dp
             in_correct_cell = in_correct_cell .and. abs(frac) < threshold
          end if
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i1 == p%icell%i1
       end if

       if(p%on_wall_id%w2 == -1) then
          frac = theta / geo%w2(p%icell%i2) - 1._dp
          in_correct_cell = in_correct_cell .and. abs(frac) < threshold
       else if(p%on_wall_id%w2 == +1) then
          frac =theta / geo%w2(p%icell%i2 + 1) - 1._dp
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
    real(dp) :: r,t,p

    call random(r)
    call random(t)
    call random(p)
    r = (r * (geo%w1(icell%i1+1)**3. - geo%w1(icell%i1)**3.) + geo%w1(icell%i1)**3.) ** (1./3.)
    t = acos(t * (geo%wcost(icell%i2+1) - geo%wcost(icell%i2)) + geo%wcost(icell%i2))
    p = p * (geo%w3(icell%i3+1) - geo%w3(icell%i3)) + geo%w3(icell%i3)

    ! If the photon doesn't fall within the walls, the cell is probably very
    ! narrow along that direction, so can just choose the mid-point

    if(r <= geo%w1(icell%i1) .or. r >= geo%w1(icell%i1+1)) then
       r = 0.5_dp * (geo%w1(icell%i1) + geo%w1(icell%i1+1))
    end if

    if(t <= geo%w2(icell%i2) .or. t >= geo%w2(icell%i2+1)) then
       t = 0.5_dp * (geo%w2(icell%i2) + geo%w2(icell%i2+1))
    end if

    if(p <= geo%w3(icell%i3) .or. p >= geo%w3(icell%i3+1)) then
       p = 0.5_dp * (geo%w3(icell%i3) + geo%w3(icell%i3+1))
    end if

    pos%x = r * sin(t) * cos(p)
    pos%y = r * sin(t) * sin(p)
    pos%z = r * cos(t)

  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p
    real(dp) :: r, rcyl
    real(dp) :: d1,d2,d3,d4,d5,d6

    ! Spherical walls - point-sphere distance

    r = sqrt(p%r .dot. p%r)
    d1 = r - geo%w1(p%icell%i1)
    d2 = geo%w1(p%icell%i1+1) - r

    ! If value is within machine precision of wall position, then set to zero
    ! (this takes care of negative values too)
    if(abs(d1) < geo%ew1(p%icell%i1)) d1 = 0._dp
    if(abs(d2) < geo%ew1(p%icell%i1+1)) d2 = 0._dp

    ! Theta walls - point-line distance
    !
    ! Distance between a line with equation y=-a/b*x-c/b and point (x0,y0) is:
    !
    ! d = |a*x0+b*y0+c|/sqrt(a**2+b**2)
    !
    ! In case of theta walls, we have y = x / tan(theta) so:
    !
    ! a = -1
    ! b = tan(theta)
    ! c = 0
    !
    ! and therefore
    !
    ! d = |-x0+tan(theta)*y0| / sqrt[1 + tan^2(theta)]

    rcyl = sqrt(p%r%x*p%r%x + p%r%y*p%r%y)
    d3 = abs(-rcyl + geo%wtant(p%icell%i2) * p%r%z) / sqrt(1 + geo%wtant(p%icell%i2)**2)
    d4 = abs(-rcyl + geo%wtant(p%icell%i2+1) * p%r%z) / sqrt(1 + geo%wtant(p%icell%i2+1)**2)

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

    real(dp) :: pA,pB,pC,pC_1,pC_2,t1,t2,z1,z2
    real(dp) :: v2_xy,v2_z
    real(dp) :: rv_xy,rv_z
    real(dp) :: r2_xy,r2_z

    real(dp) :: x_i,y_i,phi_i, dphi

    ! Find the intersections with all walls:

    call reset_t()

    ! -------------------------------------------------
    ! Compute general values
    ! -------------------------------------------------

    v2_xy = p%v%x*p%v%x+p%v%y*p%v%y
    v2_z  = p%v%z*p%v%z

    rv_xy = p%r%x*p%v%x+p%r%y*p%v%y
    rv_z  = p%r%z*p%v%z

    r2_xy = p%r%x*p%r%x+p%r%y*p%r%y
    r2_z  = p%r%z*p%r%z

    ! -------------------------------------------------
    ! Spherical walls
    ! -------------------------------------------------

    ! Pre-compute some coefficients of the quadratic that don't depend on
    ! the wall, only the position and direction of the photon.
    pB=rv_xy + rv_z ; pB = pB + pB
    pC=r2_xy + r2_z

    ! Check for intersections with inner wall. We only need to check this if
    ! we are not moving radially outwards.
    if(.not.radial) then

       ! Check for intersections with inner wall.
       pC_1 = pC - geo%wr2(p%icell%i1)
       call quadratic_pascal_reduced(pB, pC_1, t1, t2)

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
          call insert_t(t1, 1, -1, geo%ew1(p%icell%i1))
          call insert_t(t2, 1, -1, geo%ew1(p%icell%i1))
       end if

    end if

    ! Check for intersections with outer wall.
    pC_2 = pC - geo%wr2(p%icell%i1 + 1)
    call quadratic_pascal_reduced(pB, pC_2, t1, t2)

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
       call insert_t(t1, 1, +1, geo%ew1(p%icell%i1 + 1))
       call insert_t(t2, 1, +1, geo%ew1(p%icell%i1 + 1))
    end if

    ! -------------------------------------------------
    ! Cone walls
    ! -------------------------------------------------

    ! Check for intersections with lower wall. We should not check for
    ! intersections with the theta=0 wall.
    if(p%icell%i2 > 1) then

       ! Check if photon is on wall, and if so, whether it is moving along
       ! wall. If so, don't check for intersections and just set iext, which
       ! is used to specify walls that the photon is on even without an
       ! intersection.
       if(p%on_wall_id%w2 == -1 .and. equal_nulp(geo%wtant(p%icell%i2), sqrt(v2_xy) / p%v%z, 10)) then

          iext%w2 = -1

       else

          ! If the wall is the midplane wall, we can just treat it as a
          ! plane, so don't need to do any fancy calculations.
          if(p%icell%i2 == geo%midplane .and. p%v%z .ne. 0) then

             ! If the photon is on the wall, then we don't add any intersections
             if(p%on_wall_id%w2 /= -1) call insert_t(-p%r%z/p%v%z,2, -1, geo%ew2(p%icell%i2))

          else

             ! Compute the full intersections. Find the coefficients of the quadratic equation
             pA = v2_xy - v2_z * geo%wtant2(p%icell%i2)
             pB = rv_xy - rv_z * geo%wtant2(p%icell%i2) ; pB = pB + pB
             pC = r2_xy - r2_z * geo%wtant2(p%icell%i2)

             if(abs(pA) .gt. 0._dp) then  ! Solve ax^2 + bx + c = 0

                ! Compute the full quadratic solutions
                call quadratic(pA, pB, pC, t1, t2)

                ! Check if the solutions are on the right side of the
                ! mid-plane. The equation of the cone (z^2 = x^2 + y^2) does
                ! not differentiate between a positive and a negative value,
                ! so we have to check against the angle that was originally
                ! specified for the wall.
                z1 = p%r%z + p%v%z * t1
                if(z1 > 0._dp .neqv. geo%wtant(p%icell%i2) > 0._dp) t1 = huge(1._dp)
                z2 = p%r%z + p%v%z * t2
                if(z2 > 0._dp .neqv. geo%wtant(p%icell%i2) > 0._dp) t2 = huge(1._dp)

                ! If we are on the wall, then we should discard the
                ! intersection with the smallest absolute value as this will
                ! be the wall we are on. Otherwise, we should include both
                ! values.
                if(p%on_wall_id%w2 == -1) then
                   if(abs(t1) < abs(t2)) then
                      call insert_t(t2, 2, -1, geo%ew2(p%icell%i2))
                   else
                      call insert_t(t1, 2, -1, geo%ew2(p%icell%i2))
                   end if
                else
                   call insert_t(t1, 2, -1,geo%ew2(p%icell%i2))
                   call insert_t(t2, 2, -1,geo%ew2(p%icell%i2))
                end if

             else if(abs(pB) .gt. 0._dp) then  ! Solve bx + c = 0

                ! If the photon is on the wall, then we don't add any intersections
                if(p%on_wall_id%w2 /= -1) call insert_t(-pC / pB, 2, -1, geo%ew2(p%icell%i2))

             end if

          end if

       end if

    end if

    ! Check for intersections with upper wall. We should not check for
    ! intersections with the theta=pi wall.
    if(p%icell%i2 < geo%n2) then  ! originally i + 1 < n + 1, but simplify for performance

       ! Check if photon is on wall, and if so, whether it is moving along
       ! wall. If so, don't check for intersections and just set iext, which
       ! is used to specify walls that the photon is on even without an
       ! intersection.
       if(p%on_wall_id%w2 == +1 .and. equal_nulp(geo%wtant(p%icell%i2 + 1), sqrt(v2_xy) / p%v%z, 10)) then

          iext%w2 = +1

       else

          ! If the wall is the midplane wall, we can just treat it as a
          ! plane, so don't need to do any fancy calculations.
          if(p%icell%i2 + 1 == geo%midplane .and. p%v%z .ne. 0) then

             ! If the photon is on the wall, then we don't add any intersections
             if(p%on_wall_id%w2 /= +1) call insert_t(-p%r%z / p%v%z, 2, +1, geo%ew2(p%icell%i2 + 1))

          else

             ! Compute the full intersections. Find the coefficients of the quadratic equation
             pA=v2_xy - v2_z * geo%wtant2(p%icell%i2 + 1)
             pB=rv_xy - rv_z * geo%wtant2(p%icell%i2 + 1) ; pB = pB + pB
             pC=r2_xy - r2_z * geo%wtant2(p%icell%i2 + 1)

             if(abs(pA) .gt. 0._dp) then  ! Solve ax^2 + bx + c = 0

                ! Compute the full quadratic solutions
                call quadratic(pA, pB, pC, t1, t2)

                ! Check if the solutions are on the right side of the
                ! mid-plane. The equation of the cone (z^2 = x^2 + y^2) does
                ! not differentiate between a positive and a negative value,
                ! so we have to check against the angle that was originally
                ! specified for the wall.
                z1=p%r%z + p%v%z * t1
                if(z1 > 0._dp .neqv. geo%wtant(p%icell%i2 + 1) > 0._dp) t1 = huge(1._dp)
                z2=p%r%z + p%v%z * t2
                if(z2 > 0._dp .neqv. geo%wtant(p%icell%i2 + 1) > 0._dp) t2 = huge(1._dp)

                ! If we are on the wall, then we should discard the
                ! intersection with the smallest absolute value as this will
                ! be the wall we are on. Otherwise, we should include both
                ! values.
                if(p%on_wall_id%w2 == +1) then
                   if(abs(t1) < abs(t2)) then
                      call insert_t(t2, 2, +1, geo%ew2(p%icell%i2 + 1))
                   else
                      call insert_t(t1, 2, +1, geo%ew2(p%icell%i2 + 1))
                   end if
                else
                   call insert_t(t1, 2, +1,geo%ew2(p%icell%i2 + 1))
                   call insert_t(t2, 2, +1,geo%ew2(p%icell%i2 + 1))
                end if

             else if(abs(pB) .gt. 0._dp) then  ! Solve bx + c = 0

                ! If the photon is on the wall, then we don't add any intersections
                if(p%on_wall_id%w2 /= +1) call insert_t(-pC / pB, 2, +1, geo%ew2(p%icell%i2 + 1))

             end if

          end if

       end if

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
