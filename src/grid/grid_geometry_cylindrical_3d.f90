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
  public :: setup_grid_geometry

  logical :: debug = .false.

  real(dp) :: tn1,tp1,tp2
  integer :: in1,ip1,ip2

  public :: escaped
  interface escaped
     module procedure escaped_photon
     module procedure escaped_cell
  end interface escaped

  type(grid_geometry_desc),public,target :: geo

contains

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

  end subroutine setup_grid_geometry

  subroutine grid_geometry_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_geometry_debug

  type(grid_cell) function find_cell(p) result(icell)
    implicit none
    type(photon),intent(in) :: p
    real(dp) :: r_squared,phi
    integer :: i1, i2, i3
    if(debug) write(*,'(" [debug] find_cell")')
    r_squared = p%r%x*p%r%x+p%r%y*p%r%y
    phi = atan2(p%r%y,p%r%x)
    if(phi < 0._dp) phi = phi + twopi
    i1 = locate(geo%wr2,r_squared)
    i2 = locate(geo%w2,p%r%z)
    i3 = locate(geo%w3,phi)
    if(i1<1.or.i1>geo%n1) then
       call warn("find_cell","photon not in cell (in r direction)")
       icell = invalid_cell
       return
    end if
    if(i2<1.or.i2>geo%n2) then
       call warn("find_cell","photon not in cell (in z direction)")
       icell = invalid_cell
       return
    end if
    if(i3<1.or.i3>geo%n3) then
       call warn("find_cell","photon not in cell (in phi direction)")
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
    if(cell%i1<1) return
    if(cell%i1>geo%n1) return
    if(cell%i2<1) return
    if(cell%i2>geo%n2) return
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
       if(i3==0) i3 = geo%n3
    case(6)
       i3 = i3 + 1
       if(i3==geo%n3+1) i3 = 1
    end select
    next_cell = new_grid_cell(i1, i2, i3, geo)
  end function next_cell

  logical function in_correct_cell(p)
    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: rad,phi,frac,dphi
    icell_actual = find_cell(p)
    if(p%on_wall) then
       rad = sqrt(p%r%x*p%r%x+p%r%y*p%r%y)
       phi = atan2(p%r%y,p%r%x)
       if(phi < 0._dp) phi = phi + twopi
       select case(p%on_wall_id)
       case(1)
          frac = sqrt(rad / geo%w1(p%icell%i1)) - 1._dp
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(2)
          frac = sqrt(rad / geo%w1(p%icell%i1+1)) - 1._dp
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(3)
          frac = (p%r%z - geo%w2(p%icell%i2)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(4)
          frac = (p%r%z - geo%w2(p%icell%i2+1)) / (geo%w2(p%icell%i2+1) - geo%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(5)
          dphi = phi - geo%w3(p%icell%i3)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
          frac = dphi / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i2 == p%icell%i2
       case(6)
          dphi = phi - geo%w3(p%icell%i3+1)
          if(dphi > pi) dphi = dphi - twopi
          if(dphi < -pi) dphi = dphi + twopi
          frac = dphi / (geo%w3(p%icell%i3+1) - geo%w3(p%icell%i3))
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

  subroutine find_wall(p,radial,tmin,id_min)

    implicit none

    type(photon), intent(inout) :: p
    ! Position and direction

    logical,intent(in) :: radial

    integer,intent(out) :: id_min
    ! ID of next wall

    real(dp),intent(out)  :: tmin
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

    pB=rv_xy / v2_xy ; pB = pB + pB
    pC=r2_xy / v2_xy

    ! Inner wall
    pC_1 = pC - geo%wr2(p%icell%i1) / v2_xy
    call quadratic_pascal_reduced(pB,pC_1,t1,t2)
    call insert_t(t1,1)
    call insert_t(t2,1)

    ! Outer wall
    pC_2 = pC - geo%wr2(p%icell%i1+1) / v2_xy
    call quadratic_pascal_reduced(pB,pC_2,t1,t2)
    call insert_t(t1,2)
    call insert_t(t2,2)

    ! -------------------------------------------------
    ! z walls
    ! -------------------------------------------------

    t1 = ( geo%w2(p%icell%i2)   - p%r%z ) / p%v%z
    call insert_t(t1,3)
    t2 = ( geo%w2(p%icell%i2+1) - p%r%z ) / p%v%z
    call insert_t(t2,4)

    ! -------------------------------------------------
    ! phi walls
    ! -------------------------------------------------

    if(geo%n_dim == 3) then

       ! Find intersection with lower phi wall
       t1 = - ( geo%wtanp(p%icell%i3) * p%r%x - p%r%y ) / ( geo%wtanp(p%icell%i3) * p%v%x - p%v%y )

       ! Find position of intersection in x,y
       x_i = p%r%x + p%v%x * t1
       y_i = p%r%y + p%v%y * t1
       phi_i = atan2(y_i, x_i)
       dphi = abs(phi_i - geo%w3(p%icell%i3))
       if(dphi > pi) dphi = abs(dphi - twopi)
       if(dphi < 0.5 * pi) call insert_t(t1, 5)

       ! Find intersection with upper phi wall
       t2 = - ( geo%wtanp(p%icell%i3+1) * p%r%x - p%r%y ) / ( geo%wtanp(p%icell%i3+1) * p%v%x - p%v%y )

       ! Find position of intersection in x,y
       x_i = p%r%x + p%v%x * t2
       y_i = p%r%y + p%v%y * t2
       phi_i = atan2(y_i, x_i)
       dphi = abs(phi_i - geo%w3(p%icell%i3+1))
       if(dphi > pi) dphi = abs(dphi - twopi)
       if(dphi < 0.5 * pi) call insert_t(t2, 6)

    end if

    call find_next_wall(p%on_wall,tmin,id_min)

  end subroutine find_wall

  subroutine reset_t()
    implicit none
    tn1 = -huge(tn1)
    tp1 = +huge(tp1)
    tp2 = +huge(tp2)
    in1 = 0
    ip1 = 0
    ip2 = 0
  end subroutine reset_t

  subroutine insert_t(t,i)
    implicit none
    real(dp),intent(in)    :: t
    integer,intent(in)    :: i
    if(debug) print *,'[debug] inserting t,i=',t,i
    if(t < 0.) then
       if(t > tn1) then
          tn1 = t
          in1 = i
       end if
    else
       if(t < tp1) then
          tp2 = tp1
          ip2 = ip1
          tp1 = t
          ip1 = i
       else if(t < tp2) then
          tp2 = t
          ip2 = i
       end if
    end if
  end subroutine insert_t

  subroutine find_next_wall(on_wall,t,i)
    implicit none
    logical(1),intent(in)  :: on_wall
    real(dp),intent(out)    :: t
    integer,intent(out)    :: i
    if(on_wall) then
       if(abs(tn1) < abs(tp1)) then
          t = tp1
          i = ip1
       else
          t = tp2
          i = ip2
       end if
    else
       t = tp1
       i = ip1
    end if
    if(debug) print *,'[debug] selecting t,i=',t,i
  end subroutine find_next_wall

end module grid_geometry_specific
