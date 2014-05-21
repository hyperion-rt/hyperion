module grid_geometry_specific

  use core_lib
  use mpi_core
  use mpi_hdf5_io
  use type_photon
  use type_grid_cell
  use type_grid
  use grid_io
  use counters
  use kdtree2_module

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

  integer :: n_filled = 0

  ! Declare here otherwise pointer becomes corrupt
  real(dp),allocatable, target :: points(:,:)

contains

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    select case(idir)
    case(1)
       cell_width = 0.  ! TODO
    case(2)
       cell_width = 0.  ! TODO
    case(3)
       cell_width = 0.  ! TODO
    end select
  end function cell_width

  real(dp) function cell_area(cell, iface)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: iface
    select case(iface)
    case(1,2)
       cell_area = 0.  ! TODO
    case(3,4)
       cell_area = 0.  ! TODO
    case(5,6)
       cell_area = 0.  ! TODO
    end select
  end function cell_area

  ! Standard Routines

  subroutine setup_grid_geometry(group)

    ! Read a voronoi mesh from an HDF5 file

    implicit none

    integer(hid_t),intent(in) :: group

    integer :: ic, iv
    real(dp), allocatable :: x(:), y(:), z(:)
    integer, allocatable :: neighbors(:,:)
    integer :: n_neighbors
    real(dp),allocatable :: bmin(:,:), bmax(:,:)
    type(cell),pointer :: c

    ! Read geometry file
    call mp_read_keyword(group, '.', "geometry", geo%id)
    call mp_read_keyword(group, '.', "grid_type", geo%type)

    ! Read in list of cells
    call mp_table_read_column_auto(group, 'cells', 'coordinates', points)
    call mp_table_read_column_auto(group, 'cells', 'neighbours', neighbors)
    call mp_table_read_column_auto(group, 'cells', 'bmin', bmin)
    call mp_table_read_column_auto(group, 'cells', 'bmax', bmax)

    allocate(x(size(points,2)))
    allocate(y(size(points,2)))
    allocate(z(size(points,2)))

    x = points(1,:)
    y = points(2,:)
    z = points(3,:)

    ! Find number of cells
    geo%n_cells = size(x)
    geo%n_masked = geo%n_cells

    ! Allocate cells
    allocate(geo%cells(geo%n_cells))

    ! Assigned refined values to all cells
    do ic=1,geo%n_cells
       geo%cells(ic)%r%x = x(ic)
       geo%cells(ic)%r%y = y(ic)
       geo%cells(ic)%r%z = z(ic)
       n_neighbors = count(neighbors(:, ic) > -1)
       allocate(geo%cells(ic)%neighbors(n_neighbors))
       geo%cells(ic)%neighbors(1:n_neighbors) = neighbors(1:n_neighbors, ic) + 1
    end do

    ! Construct KDTree
    geo%tree => kdtree2_create(points, rearrange=.true., sort=.true.)

    ! Read parameters for top-level cell
    call mp_read_keyword(group, '.', 'xmin', geo%xmin)
    call mp_read_keyword(group, '.', 'xmax', geo%xmax)
    call mp_read_keyword(group, '.', 'ymin', geo%ymin)
    call mp_read_keyword(group, '.', 'ymax', geo%ymax)
    call mp_read_keyword(group, '.', 'zmin', geo%zmin)
    call mp_read_keyword(group, '.', 'zmax', geo%zmax)

    call mp_table_read_column_auto(group, 'cells', 'volume', geo%volume)

    geo%masked = .true.
    allocate(geo%mask(geo%n_cells))
    geo%mask = geo%volume > 0._dp
    geo%n_masked = count(geo%mask)
    allocate(geo%mask_map(geo%n_masked))
    iv = 0
    do ic=1,geo%n_cells
       if(geo%mask(ic)) then
          iv = iv + 1
          geo%mask_map(iv) = ic
       end if
    end do

    where(geo%volume < 0._dp) geo%volume = 0._dp

    geo%n_dim = 3

    ! Determine rough bounding box. Not correct, but will do the trick for now.
    geo%cells%xmin = bmin(1,:)
    geo%cells%xmax = bmax(1,:)
    geo%cells%ymin = bmin(2,:)
    geo%cells%ymax = bmax(2,:)
    geo%cells%zmin = bmin(3,:)
    geo%cells%zmax = bmax(3,:)

  end subroutine setup_grid_geometry

  subroutine grid_geometry_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_geometry_debug

  type(grid_cell) function find_cell(p) result(icell)
    implicit none
    type(photon),intent(in) :: p
    integer :: ic
    type(kdtree2_result) :: results(1)
    real(dp) :: point(3)
    if(debug) write(*,'(" [debug] find_cell")')
    ! TODO: nearest-neighbor algorithm - potentially going to be the
    ! bottleneck.
    point = [p%r%x, p%r%y, p%r%z]
    call kdtree2_n_nearest(geo%tree, point, 1, results)
    ic = results(1)%idx
    icell = new_grid_cell(ic, geo)
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
    if(cell%ic == geo%n_cells + 1) return
    escaped_cell = .false.
  end function escaped_cell

  recursive type(grid_cell) function next_cell_int(cell, direction, intersection) result(c)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    type(grid_cell) :: parent_cell
    integer :: subcell_id
    ! wall ID = cell ID
    c = new_grid_cell(direction, geo)
  end function next_cell_int

  type(grid_cell) function next_cell_wall_id(cell, direction, intersection) result(c)
    implicit none
    type(grid_cell),intent(in) :: cell
    type(wall_id),intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    c = new_grid_cell(direction%w1, geo)
  end function next_cell_wall_id

  logical function in_correct_cell(p)
    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: curr
    type(grid_cell) :: icell_actual
    real(dp) :: frac,frac1, frac2, frac3
    real(dp) :: point(3)
    type(kdtree2_result) :: results(2)
    ! TODO - only check the second result if it's close to the first
    point = [p%r%x, p%r%y, p%r%z]
    call kdtree2_n_nearest(geo%tree, point, 2, results)
    in_correct_cell = p%icell%ic == results(1)%idx .or. p%icell%ic == results(2)%idx
  end function in_correct_cell

  subroutine random_position_cell(icell,pos)
    implicit none
    type(grid_cell),intent(in) :: icell
    type(vector3d_dp), intent(out) :: pos
    type(grid_cell) :: icell_actual

    type(kdtree2_result) :: results(1)
    real(dp) :: point(3)
    integer :: i

    do i=1,1000000

        ! Sample in bounding box
        call random_uni(pos%x, geo%cells(icell%ic)%xmin, geo%cells(icell%ic)%xmax)
        call random_uni(pos%y, geo%cells(icell%ic)%ymin, geo%cells(icell%ic)%ymax)
        call random_uni(pos%z, geo%cells(icell%ic)%zmin, geo%cells(icell%ic)%zmax)

        ! Check if in right cell
        point = [pos%x, pos%y, pos%z]
        call kdtree2_n_nearest(geo%tree, point, 1, results)

        ! Break if correct
        if(results(1)%idx == icell%ic) exit

    end do

    if(i==1000001) call error("random_position_cell", "too many samples")

  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)
    implicit none
    type(photon),intent(in) :: p
    real(dp) :: d1,d2,d3,d4,d5,d6
    ! Is this the second nearest neighbor?
  end function distance_to_closest_wall

  subroutine find_wall(p,radial,tmin,id_min)

    implicit none

    type(photon), intent(inout) :: p
    ! Position and direction

    logical,intent(in) :: radial

    type(wall_id),intent(out) :: id_min
    ! ID of next wall

    real(dp),intent(out)  :: tmin
    ! tmin to nearest wall

    type(cell),pointer :: icell, ocell


    integer :: i, n_neighbors
    type(vector3d_dp) :: r0, ri, rmid, n
    real(dp) :: t, tx, ty, tz


    icell => geo%cells(p%icell%ic)

    n_neighbors = size(icell%neighbors)

    ! Now we cycle through all the neighbors and find for each the
    ! intersection with the plane bisecting the points.

    ! First find intersection with outer boundary

    call reset_t()

    if(p%v%x > 0._dp) then
        tx = ( geo%xmax - p%r%x ) / p%v%x
    else if(p%v%x < 0._dp) then
        tx = ( geo%xmin - p%r%x) / p%v%x
    else
        tx = huge(1._dp)
    end if
    call insert_t(tx, 1, geo%n_cells + 1, 0._dp)

    if(p%v%y > 0.) then
        ty = ( geo%ymax - p%r%y ) / p%v%y
    else if(p%v%y < 0._dp) then
        ty = ( geo%ymin - p%r%y ) / p%v%y
    else
        ty = huge(1._dp)
    end if
    call insert_t(ty, 1, geo%n_cells + 1, 0._dp)

    if(p%v%z > 0.) then
        tz = ( geo%zmax - p%r%z ) / p%v%z
    else if(p%v%z < 0._dp) then
        tz = ( geo%zmin - p%r%z ) / p%v%z
    else
        tz = huge(1._dp)
    end if
    call insert_t(tz, 1, geo%n_cells + 1, 0._dp)

    do i = 1,n_neighbors

        if(icell%neighbors(i) == -p%on_wall_id%w2) cycle

        ! Declare pointer to other cell
        ocell => geo%cells(icell%neighbors(i))

        n = ocell%r - icell%r
        rmid = 0.5_dp * (ocell%r + icell%r)

        t = (n.dot.(rmid - p%r)) / (n.dot.p%v)

        call insert_t(t, 1, icell%neighbors(i), 0._dp)

    end do

    call find_next_wall(tmin,id_min)

    ! use w2 to store previous cell
    id_min%w2 = p%icell%ic

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
    i = imin
    if(debug) print *,'[debug] selecting t,i=',t,i
  end subroutine find_next_wall

end module grid_geometry_specific
