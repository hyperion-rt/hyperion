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

  ! Oct cell organization:
  ! subcell 1: (0,0,0)
  ! subcell 2: (1,0,0)
  ! subcell 3: (0,1,0)
  ! subcell 4: (1,1,0)
  ! subcell 5: (0,0,1)
  ! subcell 6: (1,0,1)
  ! subcell 7: (0,1,1)
  ! subcell 8: (1,1,1)

  ! The following array gives the opposite cell for a given subcell and
  ! wall id. The call signature is (subcell_id, wall_id)
  integer,parameter :: opposite_cell(8,6) = reshape((/0,1,0,3,0,5,0,7,&
       &                                              2,0,4,0,6,0,8,0,&
       &                                              0,0,1,2,0,0,5,6,&
       &                                              3,4,0,0,7,8,0,0,&
       &                                              0,0,0,0,1,2,3,4,&
       &                                              5,6,7,8,0,0,0,0&
       &                                              /),(/8,6/))

  public :: setup_grid_geometry
  type(grid_geometry_desc),public,target :: geo

  integer :: n_filled = 0

contains

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    select case(idir)
    case(1)
       cell_width = geo%cells(cell%ic)%dx * 2._dp
    case(2)
       cell_width = geo%cells(cell%ic)%dy * 2._dp
    case(3)
       cell_width = geo%cells(cell%ic)%dz * 2._dp
    end select
  end function cell_width

  real(dp) function cell_area(cell, iface)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: iface
    select case(iface)
    case(1,2)
       cell_area = geo%cells(cell%ic)%dy * geo%cells(cell%ic)%dz * 4._dp
    case(3,4)
       cell_area = geo%cells(cell%ic)%dz * geo%cells(cell%ic)%dx * 4._dp
    case(5,6)
       cell_area = geo%cells(cell%ic)%dx * geo%cells(cell%ic)%dy * 4._dp
    end select
  end function cell_area

  ! Octree Helper Routines

  integer function subcell_id(cell_id, r)
    ! Given a cell and a position falling in the cell, find the sub-cell
    ! that the position would fall into
    implicit none
    integer,intent(in) :: cell_id
    type(vector3d_dp),intent(in) :: r
    if(r%x < geo%cells(cell_id)%x) then
       if(r%y < geo%cells(cell_id)%y) then
          if(r%z < geo%cells(cell_id)%z) then
             subcell_id = 1
          else
             subcell_id = 5
          end if
       else
          if(r%z < geo%cells(cell_id)%z) then
             subcell_id = 3
          else
             subcell_id = 7
          end if
       end if
    else
       if(r%y < geo%cells(cell_id)%y) then
          if(r%z < geo%cells(cell_id)%z) then
             subcell_id = 2
          else
             subcell_id = 6
          end if
       else
          if(r%z < geo%cells(cell_id)%z) then
             subcell_id = 4
          else
             subcell_id = 8
          end if
       end if
    end if
  end function subcell_id

  recursive integer function locate_cell(r, cell_id) result(final_cell_id)
    implicit none
    type(vector3d_dp), intent(in) :: r
    integer,intent(in) :: cell_id
    integer :: new_cell_id
    if(geo%cells(cell_id)%refined) then
       new_cell_id = geo%cells(cell_id)%children(subcell_id(cell_id, r))
       final_cell_id = locate_cell(r, new_cell_id)
    else
       final_cell_id = cell_id
    end if
  end function locate_cell

  recursive subroutine octree_setup_indiv(parent_id)

    ! Set up a cell that has already had the refinement flag set

    implicit none

    integer,intent(in) :: parent_id
    integer :: sx, sy, sz
    integer :: ic, child_id

    allocate(geo%cells(parent_id)%children(8))

    do ic=1,8

       n_filled = n_filled + 1
       child_id = n_filled
       geo%cells(parent_id)%children(ic) = child_id

       sx = 1 ; sy = 1 ; sz = 1
       if(mod(ic-1,2)==0) sx = -sx
       if(mod(ic-1-mod(ic-1,2),4)/2==0) sy = -sy
       if(mod(ic-1-mod(ic-1,2)-mod(ic-1-mod(ic-1,2),4),8)/4==0) sz = -sz

       geo%cells(child_id)%x = geo%cells(parent_id)%x + sx * geo%cells(parent_id)%dx/2._8
       geo%cells(child_id)%y = geo%cells(parent_id)%y + sy * geo%cells(parent_id)%dy/2._8
       geo%cells(child_id)%z = geo%cells(parent_id)%z + sz * geo%cells(parent_id)%dz/2._8
       geo%cells(child_id)%dx = geo%cells(parent_id)%dx / 2._dp
       geo%cells(child_id)%dy = geo%cells(parent_id)%dy / 2._dp
       geo%cells(child_id)%dz = geo%cells(parent_id)%dz / 2._dp

       geo%cells(child_id)%parent = parent_id
       geo%cells(child_id)%parent_subcell = ic

       if(geo%cells(child_id)%refined) call octree_setup_indiv(child_id)

    end do

  end subroutine octree_setup_indiv

  ! Standard Routines

  subroutine setup_grid_geometry(group)

    ! Read an octree from an HDF5 file

    implicit none

    integer(hid_t),intent(in) :: group

    real(dp) :: x, y, z, dx, dy, dz
    integer :: ic, iv

    ! FIX - issue with reading 64-bit 1D column with 32-bit variables
    integer, allocatable :: refined(:)

    ! Read geometry file
    call mp_read_keyword(group, '.', "geometry", geo%id)
    call mp_read_keyword(group, '.', "grid_type", geo%type)

    ! Read in refinement list
    call mp_table_read_column_auto(group, 'cells', 'refined', refined)

    ! Find number of cells
    geo%n_cells = size(refined)

    ! Figure out which cells are valid, and construct an index to map valid
    ! cell IDs to original IDs
    geo%masked = .true.
    allocate(geo%mask(geo%n_cells))
    geo%mask = refined == 0
    geo%n_masked = count(geo%mask)
    allocate(geo%mask_map(geo%n_masked))
    iv = 0
    do ic=1,geo%n_cells
       if(geo%mask(ic)) then
          iv = iv + 1
          geo%mask_map(iv) = ic
       end if
    end do

    ! Allocate cells
    allocate(geo%cells(geo%n_cells))

    ! Assigned refined values to all cells
    do ic=1,geo%n_cells
       geo%cells(ic)%refined = refined(ic) == 1
    end do

    ! Read parameters for top-level cell
    call mp_read_keyword(group, '.', 'x', geo%cells(1)%x)
    call mp_read_keyword(group, '.', 'y', geo%cells(1)%y)
    call mp_read_keyword(group, '.', 'z', geo%cells(1)%z)
    call mp_read_keyword(group, '.', 'dx', geo%cells(1)%dx)
    call mp_read_keyword(group, '.', 'dy', geo%cells(1)%dy)
    call mp_read_keyword(group, '.', 'dz', geo%cells(1)%dz)

    n_filled=1

    ! Recursively set up cells
    if(geo%cells(1)%refined) call octree_setup_indiv(1)

    allocate(geo%volume(geo%n_cells))

    ! Compute cell volumes
    do ic=1,geo%n_cells
       geo%volume(ic) = geo%cells(ic)%dx * geo%cells(ic)%dy * geo%cells(ic)%dz * 8._dp
    end do

    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')

    geo%n_dim = 3

    geo%xmin = geo%cells(1)%x - geo%cells(1)%dx
    geo%xmax = geo%cells(1)%x + geo%cells(1)%dx
    geo%ymin = geo%cells(1)%y - geo%cells(1)%dy
    geo%ymax = geo%cells(1)%y + geo%cells(1)%dy
    geo%zmin = geo%cells(1)%z - geo%cells(1)%dz
    geo%zmax = geo%cells(1)%z + geo%cells(1)%dz

    ! Set precision
    geo%eps = spacing(max(geo%cells(1)%dx, geo%cells(1)%dy, geo%cells(1)%dz)) * 3._dp

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
    if(debug) write(*,'(" [debug] find_cell")')
    if(p%r%x<geo%xmin.or.p%r%y>geo%xmax) then
       call warn("find_cell","photon not in grid (in x direction)")
       icell = invalid_cell
       return
    end if
    if(p%r%y<geo%ymin.or.p%r%y>geo%ymax) then
       call warn("find_cell","photon not in grid (in y direction)")
       icell = invalid_cell
       return
    end if
    if(p%r%z<geo%zmin.or.p%r%y>geo%zmax) then
       call warn("find_cell","photon not in grid (in z direction)")
       icell = invalid_cell
       return
    end if
    ic = locate_cell(p%r, 1)
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
    if(cell%ic==1) then
       c = new_grid_cell(geo%n_cells+1, geo)
       return
    end if
    subcell_id = opposite_cell(geo%cells(cell%ic)%parent_subcell, direction)
    parent_cell = new_grid_cell(geo%cells(cell%ic)%parent, geo)
    if(subcell_id > 0) then
       c = new_grid_cell(locate_cell(intersection, geo%cells(parent_cell%ic)%children(subcell_id)),geo)
    else
       c = next_cell_int(parent_cell, direction, intersection)
       return
    end if
  end function next_cell_int

  type(grid_cell) function next_cell_wall_id(cell, direction, intersection)
    implicit none
    type(grid_cell),intent(in) :: cell
    type(wall_id),intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    if(direction%w1 == -1) then
       next_cell_wall_id = next_cell_int(cell, 1, intersection)
    else if(direction%w1 == +1) then
       next_cell_wall_id = next_cell_int(cell, 2, intersection)
    else if(direction%w2 == -1) then
       next_cell_wall_id = next_cell_int(cell, 3, intersection)
    else if(direction%w2 == +1) then
       next_cell_wall_id = next_cell_int(cell, 4, intersection)
    else if(direction%w3 == -1) then
       next_cell_wall_id = next_cell_int(cell, 5, intersection)
    else if(direction%w3 == +1) then
       next_cell_wall_id = next_cell_int(cell, 6, intersection)
    end if
  end function next_cell_wall_id

  logical function in_correct_cell(p)
    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: curr
    type(grid_cell) :: icell_actual
    real(dp) :: frac,frac1, frac2, frac3
    icell_actual = find_cell(p)
    if(p%on_wall) then
       frac1 = abs(p%r%x - geo%cells(p%icell%ic)%x) / geo%cells(p%icell%ic)%dx
       frac2 = abs(p%r%y - geo%cells(p%icell%ic)%y) / geo%cells(p%icell%ic)%dy
       frac3 = abs(p%r%z - geo%cells(p%icell%ic)%z) / geo%cells(p%icell%ic)%dz
       if(abs(p%on_wall_id%w1) == 1) then
          frac = frac1 - 1._dp
          in_correct_cell = frac2 < 1._dp .and. frac3 < 1._dp
       end if
       if(abs(p%on_wall_id%w2) == 1) then
          frac = frac2 - 1._dp
          in_correct_cell = frac1 < 1._dp .and. frac3 < 1._dp
       end if
       if(abs(p%on_wall_id%w3) == 1) then
          frac = frac3 - 1._dp
          in_correct_cell = frac1 < 1._dp .and. frac2 < 1._dp
       end if
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
    pos%x = (2._dp*x - 1._dp) * geo%cells(icell%ic)%dx + geo%cells(icell%ic)%x
    pos%y = (2._dp*y - 1._dp) * geo%cells(icell%ic)%dy + geo%cells(icell%ic)%y
    pos%z = (2._dp*z - 1._dp) * geo%cells(icell%ic)%dz + geo%cells(icell%ic)%z
  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p

    real(dp) :: d1,d2,d3,d4,d5,d6

    d1 = p%r%x - geo%cells(p%icell%ic)%x - geo%cells(p%icell%ic)%dx
    d2 = geo%cells(p%icell%ic)%x + geo%cells(p%icell%ic)%dx - p%r%x
    d3 = p%r%y - geo%cells(p%icell%ic)%y - geo%cells(p%icell%ic)%dy
    d4 = geo%cells(p%icell%ic)%y + geo%cells(p%icell%ic)%dy - p%r%y
    d5 = p%r%z - geo%cells(p%icell%ic)%z - geo%cells(p%icell%ic)%dz
    d6 = geo%cells(p%icell%ic)%z + geo%cells(p%icell%ic)%dz - p%r%z

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

    type(wall_id),intent(out) :: id_min
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
       tx = ( geo%cells(p%icell%ic)%x + geo%cells(p%icell%ic)%dx - p%r%x ) / p%v%x
    else if(p%v%x < 0._dp) then
       tx = ( geo%cells(p%icell%ic)%x - geo%cells(p%icell%ic)%dx - p%r%x ) / p%v%x
    else
       tx = huge(1._dp)
    end if

    if(pos_vy) then
       ty = ( geo%cells(p%icell%ic)%y + geo%cells(p%icell%ic)%dy - p%r%y ) / p%v%y
    else if(p%v%y < 0._dp) then
       ty = ( geo%cells(p%icell%ic)%y - geo%cells(p%icell%ic)%dy - p%r%y ) / p%v%y
    else
       ty = huge(1._dp)
    end if

    if(pos_vz) then
       tz = ( geo%cells(p%icell%ic)%z + geo%cells(p%icell%ic)%dz - p%r%z ) / p%v%z
    else if(p%v%z < 0._dp) then
       tz = ( geo%cells(p%icell%ic)%z - geo%cells(p%icell%ic)%dz - p%r%z ) / p%v%z
    else
       tz = huge(1._dp)
    end if

    ! Find the closest of the three walls. The following effectively
    ! finds the minimum of three values. A lot of code for such a
    ! small thing, but this runs much much faster than using a built
    ! in min function or any kind of loop. Each iteraction comprises only
    ! three if statements and two pointer assignements.

    if(tx.lt.tz) then
       if(tx.lt.ty) then
          if(pos_vx) then
             id_min%w1 = +1
          else
             id_min%w1 = -1
          end if
          tmin = tx
       else
          if(pos_vy) then
             id_min%w2 = +1
          else
             id_min%w2 = -1
          end if
          tmin = ty
       end if
    else
       if(tz.lt.ty) then
          if(pos_vz) then
             id_min%w3 = +1
          else
             id_min%w3 = -1
          end if
          tmin = tz
       else
          if(pos_vy) then
             id_min%w2 = +1
          else
             id_min%w2 = -1
          end if
          tmin = ty
       end if
    end if

    if(tmin < 0._dp) then
        if(tmin > -10 * geo%eps) then
            ! TODO: there may be a better way to avoid this kind of situation
            tmin = 0._dp
        else
            call warn("find_wall", "an error occured when searching for the nearest wall (negative t)")
            id_min = no_wall
        end if
    end if

  end subroutine find_wall

end module grid_geometry_specific
