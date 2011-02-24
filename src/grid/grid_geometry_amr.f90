module grid_geometry_specific

  use core_lib
  use mpi_core
  use type_photon
  use type_grid_cell
  use type_grid
  use grid_io

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

  ! AMR Helper functions

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    cell_width = geo%levels(cell%ilevel)%fabs(cell%ifab)%width(idir)
  end function cell_width

  logical function fabs_intersect(fab1, fab2)
    ! Given two fabs, do they intersect anywhere?
    implicit none
    type(fab_desc),intent(in) :: fab1, fab2
    fabs_intersect = .false.
    if(fab1%xmax < fab2%xmin) return
    if(fab1%xmin > fab2%xmax) return
    if(fab1%ymax < fab2%ymin) return
    if(fab1%ymin > fab2%ymax) return
    if(fab1%zmax < fab2%zmin) return
    if(fab1%zmin > fab2%zmax) return
    fabs_intersect = .true.
  end function fabs_intersect

  logical function fabs_close(fab1, fab2)
    ! Given two fabs, are they within one half cell from each other
    type(fab_desc),intent(in) :: fab1, fab2
    fabs_close = .false.
    if(fab1%xmax < fab2%xmin - fab2%width(1)*0.5_dp) return
    if(fab1%xmin > fab2%xmax + fab2%width(1)*0.5_dp) return
    if(fab1%ymax < fab2%ymin - fab2%width(2)*0.5_dp) return
    if(fab1%ymin > fab2%ymax + fab2%width(2)*0.5_dp) return
    if(fab1%zmax < fab2%zmin - fab2%width(3)*0.5_dp) return
    if(fab1%zmin > fab2%zmax + fab2%width(3)*0.5_dp) return
    fabs_close = .true.
  end function fabs_close

  logical function in_fab(fab, r)
    ! Given a fab and a position r, determines if the position lies inside the fab.
    ! Could be optimized by finding distance from center (abs(x)-x0) < dx, although not clear which is faster.
    implicit none
    type(fab_desc),intent(in) :: fab
    type(vector3d_dp),intent(in) :: r
    in_fab = .false.
    if(r%x < fab%xmin) return
    if(r%x > fab%xmax) return
    if(r%y < fab%ymin) return
    if(r%y > fab%ymax) return
    if(r%z < fab%zmin) return
    if(r%z > fab%zmax) return
    in_fab = .true.
  end function in_fab

  integer function locate_fab(level, r) result(ifab)
    ! Given a level and a position r, determines if the position lies in
    ! any of the fabs, and if so returns the fab ID in the level. If the
    ! position does not lie in any fab, the function returns -1.
    implicit none
    type(level_desc),intent(in) :: level
    type(vector3d_dp),intent(in) :: r
    do ifab=1,size(level%fabs)
       if(in_fab(level%fabs(ifab), r)) return
    end do
    ifab = -1
  end function locate_fab

  subroutine read_fab(group, fab)

    ! Given and HDF5 group containing a fab, this subroutine reads and returns it

    implicit none

    integer(hid_t),intent(in) :: group
    type(fab_desc),intent(out) :: fab

    call hdf5_read_keyword(group, '.', 'n1', fab%n1)
    call hdf5_read_keyword(group, '.', 'n2', fab%n2)
    call hdf5_read_keyword(group, '.', 'n3', fab%n3)
    fab%n_cells = fab%n1 * fab%n2 * fab%n3
    fab%n_dim = 3

    allocate(fab%w1(fab%n1+1))
    call hdf5_read_keyword(group, '.', 'xmin', fab%xmin)
    call hdf5_read_keyword(group, '.', 'xmax', fab%xmax)
    call linspace(fab%xmin, fab%xmax, fab%w1)

    allocate(fab%w2(fab%n2+1))
    call hdf5_read_keyword(group, '.', 'ymin', fab%ymin)
    call hdf5_read_keyword(group, '.', 'ymax', fab%ymax)
    call linspace(fab%ymin, fab%ymax, fab%w2)

    allocate(fab%w3(fab%n3+1))
    call hdf5_read_keyword(group, '.', 'zmin', fab%zmin)
    call hdf5_read_keyword(group, '.', 'zmax', fab%zmax)
    call linspace(fab%zmin, fab%zmax, fab%w3)

    fab%width(1) = (fab%xmax-fab%xmin)/real(fab%n1, dp)
    fab%width(2) = (fab%ymax-fab%ymin)/real(fab%n2, dp)
    fab%width(3) = (fab%zmax-fab%zmin)/real(fab%n3, dp)

    fab%area(1:2) = fab%width(2) * fab%width(3)
    fab%area(3:4) = fab%width(1) * fab%width(3)
    fab%area(5:6) = fab%width(1) * fab%width(2)

    fab%volume = fab%width(1) * fab%width(2) * fab%width(3)

    allocate(fab%goto_fab(0:fab%n1+1, 0:fab%n2+1, 0:fab%n3+1))
    allocate(fab%goto_level(0:fab%n1+1, 0:fab%n2+1, 0:fab%n3+1))

    fab%goto_fab = 0
    fab%goto_level = 0

  end subroutine read_fab

  subroutine read_level(group, level)

    ! Given and HDF5 group containing a level, this subroutine reads and returns it

    integer(hid_t),intent(in) :: group
    type(level_desc),intent(out) :: level
    integer :: nfabs, ifab
    integer(hid_t) :: g_fab
    character(len=100) :: fab_name

    call hdf5_read_keyword(group, '.', 'nfabs', nfabs)

    allocate(level%fabs(nfabs))

    do ifab=1,nfabs
       write(fab_name, '("Fab ", I0)') ifab
       g_fab = hdf5_open_group(group, fab_name)
       call read_fab(g_fab, level%fabs(ifab))
       call hdf5_close_group(g_fab)
    end do

  end subroutine read_level

  ! Standard Grid Methods

  subroutine setup_grid_geometry(group)

    implicit none

    integer(hid_t),intent(in) :: group

    integer :: nlevels, ilevel, ilevel1, ilevel2
    integer(hid_t) :: g_level
    character(len=100) :: level_name

    integer :: ifab, ifab1, ifab2, icell

    integer :: start_id
    type(level_desc), pointer :: level, level1, level2
    type(fab_desc), pointer :: fab, fab1, fab2
    integer :: i1, i2, i3
    type(vector3d_dp) :: r
    real(dp) :: min_width

    ! Read geometry file
    call hdf5_read_keyword(group, '.', "geometry", geo%id)
    call hdf5_read_keyword(group, '.', "grid_type", geo%type)

    if(trim(geo%type).ne.'amr') call error("setup_grid_geometry","grid is not AMR cartesian")

    if(main_process()) write(*,'(" [setup_grid_geometry] Reading AMR cartesian grid")')

    ! Read the number of levels
    call hdf5_read_keyword(group, '.', 'nlevels', nlevels)

    ! Allocate the levels
    allocate(geo%levels(nlevels))

    ! Loop through the levels and read all the fabs
    do ilevel=1,nlevels
       write(level_name, '("Level ", I0)') ilevel
       g_level = hdf5_open_group(group, level_name)
       call read_level(g_level, geo%levels(ilevel))
       call hdf5_close_group(g_level)
    end do

    ! Find the unique ID of the first cell in each fab
    start_id = 1
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          fab%start_id = start_id
          start_id = start_id + fab%n_cells
       end do
    end do

    ! Set the total number of cells
    geo%n_cells = start_id - 1

    ! Compute volume for all cells
    allocate(geo%volume(geo%n_cells))
    min_width = huge(1._dp)
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          if(any(fab%width < min_width)) min_width = minval(fab%width)
          do icell=fab%start_id,fab%start_id+fab%n_cells-1
             geo%volume(icell) = fab%volume
             if(geo%volume(icell)==0.) stop
          end do
       end do
    end do
    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')

    ! Set number of dimensions (is this needed by PDA?)
    geo%n_dim = 3

    ! Set the step size for stepping out of fabs. Choose half the smallest width in the grid.
    geo%eps = min_width / 2._dp

    ! Find the level ID, fab ID, and internal 3D coordinates for each unique ID
    call preset_cell_id(geo)

    ! For each fab, find all cells in level below that overlap, and set flag in level below to indicate this

    do ilevel1=size(geo%levels)-1, 1, -1
       level1 => geo%levels(ilevel1)
       level2 => geo%levels(ilevel1+1)
       do ifab1=1,size(level1%fabs)
          fab1 => level1%fabs(ifab1)
          do ifab2=1,size(level2%fabs)
             fab2 => level2%fabs(ifab2)
             if(fabs_intersect(fab1, fab2)) then
                do i1=1,fab1%n1
                   r%x = 0.5_dp * (fab1%w1(i1) + fab1%w1(i1+1))
                   do i2=1,fab1%n2
                      r%y = 0.5_dp * (fab1%w2(i2) + fab1%w2(i2+1))
                      do i3=1,fab1%n3
                         r%z = 0.5_dp * (fab1%w3(i3) + fab1%w3(i3+1))
                         if(in_fab(fab2, r)) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel1+1
                         end if
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

    ! Now for all neighboring cells, find out what fab they end up in if they go one step outside the grid

    do ilevel1 = 1, size(geo%levels)
       level1 => geo%levels(ilevel1)
       do ifab1=1,size(level1%fabs)
          fab1 => level1%fabs(ifab1)

          do ilevel2 = ilevel1, 1, -1
             level2 => geo%levels(ilevel2)
             do ifab2=1,size(level2%fabs)
                fab2 => level2%fabs(ifab2)

                if(fabs_close(fab1, fab2).and.ifab1.ne.ifab2) then

                   ! xmin side
                   i1 = 0
                   r%x = fab1%xmin - fab1%width(1) * 0.5_dp
                   do i2=1,fab1%n2
                      r%y = 0.5_dp * (fab1%w2(i2) + fab1%w2(i2+1))
                      do i3=1,fab1%n3
                         r%z = 0.5_dp * (fab1%w3(i3) + fab1%w3(i3+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! xmax side
                   i1 = fab1%n1+1
                   r%x = fab1%xmax + fab1%width(1) * 0.5_dp
                   do i2=1,fab1%n2
                      r%y = 0.5_dp * (fab1%w2(i2) + fab1%w2(i2+1))
                      do i3=1,fab1%n3
                         r%z = 0.5_dp * (fab1%w3(i3) + fab1%w3(i3+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymin side
                   i2 = 0
                   r%y = fab1%ymin - fab1%width(2) * 0.5_dp
                   do i1=1,fab1%n1
                      r%x = 0.5_dp * (fab1%w1(i1) + fab1%w1(i1+1))
                      do i3=1,fab1%n3
                         r%z = 0.5_dp * (fab1%w3(i3) + fab1%w3(i3+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymax side
                   i2 = fab1%n2+1
                   r%y = fab1%ymax + fab1%width(2) * 0.5_dp
                   do i1=1,fab1%n1
                      r%x = 0.5_dp * (fab1%w1(i1) + fab1%w1(i1+1))
                      do i3=1,fab1%n3
                         r%z = 0.5_dp * (fab1%w3(i3) + fab1%w3(i3+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! zmin side
                   i3 = 0
                   r%z = fab1%zmin - fab1%width(3) * 0.5_dp
                   do i1=1,fab1%n1
                      r%x = 0.5_dp * (fab1%w1(i1) + fab1%w1(i1+1))
                      do i2=1,fab1%n2
                         r%y = 0.5_dp * (fab1%w2(i2) + fab1%w2(i2+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymax side
                   i3 = fab1%n3+1
                   r%z = fab1%zmax + fab1%width(3) * 0.5_dp
                   do i1=1,fab1%n1
                      r%x = 0.5_dp * (fab1%w1(i1) + fab1%w1(i1+1))
                      do i2=1,fab1%n2
                         r%y = 0.5_dp * (fab1%w2(i2) + fab1%w2(i2+1))
                         if(in_fab(fab2, r) .and. fab1%goto_fab(i1, i2, i3) == 0) then
                            fab1%goto_fab(i1, i2, i3) = ifab2
                            fab1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                end if

             end do
          end do

       end do
    end do

  end subroutine setup_grid_geometry

  integer function ipos2(xmin, xmax, x, nbin)
    implicit none
    real(dp),intent(in) :: xmin, xmax, x
    integer,intent(in) :: nbin
    real(dp) :: eps
    eps = (xmax-xmin) * 1.e-10
    ipos2 = ipos(xmin, xmax, x, nbin)
    if(ipos2==0 .and. abs(x - xmin) < eps) ipos2 = 1
    if(ipos2==nbin+1 .and. abs(x - xmax) < eps) ipos2 = nbin
  end function ipos2

  recursive type(grid_cell) function find_position_in_fab(r, ilevel, ifab) result(cell)
    implicit none
    type(vector3d_dp),intent(in) :: r
    integer,intent(in) :: ilevel, ifab
    integer :: i1, i2, i3, ilevel_new, ifab_new
    type(fab_desc), pointer :: fab
    fab => geo%levels(ilevel)%fabs(ifab)
    i1 = ipos2(fab%xmin, fab%xmax, r%x, fab%n1)
    i2 = ipos2(fab%ymin, fab%ymax, r%y, fab%n2)
    i3 = ipos2(fab%zmin, fab%zmax, r%z, fab%n3)
    ilevel_new = fab%goto_level(i1, i2, i3)
    ifab_new = fab%goto_fab(i1, i2, i3)
    if(ilevel_new == 0) then
       cell = new_grid_cell(i1, i2, i3, ilevel, ifab, geo)
    else
       cell = find_position_in_fab(r, ilevel_new, ifab_new)
    end if
  end function find_position_in_fab


  subroutine grid_geometry_debug(debug_flag)
    implicit none
    logical,intent(in) :: debug_flag
    debug = debug_flag
  end subroutine grid_geometry_debug

  type(grid_cell) function find_cell(p) result(icell)
    implicit none
    type(photon),intent(in) :: p
    if(debug) write(*,'(" [debug] find_cell")')    
    icell = find_cell_position(p%r)    
  end function find_cell

  type(grid_cell) function find_cell_position(r) result(icell)
    implicit none
    type(vector3d_dp),intent(in) :: r
    integer :: ifab
    if(debug) write(*,'(" [debug] find_cell")')
    ifab = locate_fab(geo%levels(1), r) ! Find fab in level 1
    icell = find_position_in_fab(r, 1, ifab) ! Refine position
  end function find_cell_position

  subroutine place_in_cell(p)
    implicit none
    type(photon),intent(inout) :: p
    p%icell = find_cell(p)
    if(p%icell == invalid_cell.or.p%icell == outside_cell) then
       call warn("place_in_cell","place_in_cell failed - killing")
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
    escaped_cell = cell == outside_cell
  end function escaped_cell

  type(grid_cell) function next_cell(cell, direction, intersection)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    type(vector3d_dp) :: r

    integer :: i1, i2, i3, ilevel, ifab
    type(fab_desc), pointer :: fab

    fab => geo%levels(cell%ilevel)%fabs(cell%ifab)

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

    if(fab%goto_level(i1, i2, i3) == 0) then
       if(i1==0.or.i1==fab%n1+1.or.i2==0.or.i2==fab%n2+1.or.i3==0.or.i3==fab%n3+1) then
          next_cell = outside_cell
       else
          next_cell = new_grid_cell(i1, i2, i3, cell%ilevel, cell%ifab, geo)
       end if
    else
       ilevel = fab%goto_level(i1, i2, i3)
       ifab = fab%goto_fab(i1, i2, i3)
       r = intersection
       select case(direction)
       case(1)
          r%x = r%x - geo%eps
       case(2)
          r%x = r%x + geo%eps
       case(3)
          r%y = r%y - geo%eps
       case(4)
          r%y = r%y + geo%eps
       case(5)
          r%z = r%z - geo%eps
       case(6)
          r%z = r%z + geo%eps
       end select
       next_cell = find_position_in_fab(r, ilevel, ifab)
    end if
  end function next_cell

  logical function in_correct_cell(p)
    ! Numerical issue here with photons that are currently on walls. Don't try and find right level and icell, assume they are correct.

    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: frac
    type(fab_desc), pointer :: fab
    fab => geo%levels(p%icell%ilevel)%fabs(p%icell%ifab)
    icell_actual = p%icell
    icell_actual%i1 = ipos(fab%xmin, fab%xmax, p%r%x, fab%n1)
    icell_actual%i2 = ipos(fab%ymin, fab%ymax, p%r%y, fab%n2)
    icell_actual%i3 = ipos(fab%zmin, fab%zmax, p%r%z, fab%n3)
    if(p%on_wall) then
       select case(p%on_wall_id)
       case(1)
          frac = (p%r%x - fab%w1(p%icell%i1)) / (fab%w1(p%icell%i1+1) - fab%w1(p%icell%i1))
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(2)
          frac = (p%r%x - fab%w1(p%icell%i1+1)) / (fab%w1(p%icell%i1+1) - fab%w1(p%icell%i1))
          in_correct_cell = icell_actual%i2 == p%icell%i2 .and. icell_actual%i3 == p%icell%i3
       case(3)
          frac = (p%r%y - fab%w2(p%icell%i2)) / (fab%w2(p%icell%i2+1) - fab%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(4)
          frac = (p%r%y - fab%w2(p%icell%i2+1)) / (fab%w2(p%icell%i2+1) - fab%w2(p%icell%i2))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i3 == p%icell%i3
       case(5)
          frac = (p%r%z - fab%w3(p%icell%i3)) / (fab%w3(p%icell%i3+1) - fab%w3(p%icell%i3))
          in_correct_cell = icell_actual%i1 == p%icell%i1 .and. icell_actual%i2 == p%icell%i2
       case(6)
          frac = (p%r%z - fab%w3(p%icell%i3+1)) / (fab%w3(p%icell%i3+1) - fab%w3(p%icell%i3))
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
    type(fab_desc), pointer :: fab
    fab => geo%levels(icell%ilevel)%fabs(icell%ifab)
    call random(x)
    call random(y)
    call random(z)
    pos%x = x * (fab%w1(icell%i1+1) - fab%w1(icell%i1)) + fab%w1(icell%i1)
    pos%y = y * (fab%w2(icell%i2+1) - fab%w2(icell%i2)) + fab%w2(icell%i2)
    pos%z = z * (fab%w3(icell%i3+1) - fab%w3(icell%i3)) + fab%w3(icell%i3)
  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p

    real(dp) :: d1,d2,d3,d4,d5,d6

    type(fab_desc), pointer :: fab

    fab => geo%levels(p%icell%ilevel)%fabs(p%icell%ifab)

    d1 = p%r%x - fab%w1(p%icell%i1)
    d2 = fab%w1(p%icell%i1+1) - p%r%x
    d3 = p%r%y - fab%w2(p%icell%i2)
    d4 = fab%w2(p%icell%i2+1) - p%r%y
    d5 = p%r%z - fab%w3(p%icell%i3)
    d6 = fab%w3(p%icell%i3+1) - p%r%z

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

    type(fab_desc), pointer :: fab

    fab => geo%levels(p%icell%ilevel)%fabs(p%icell%ifab)

    ! Can store this in photon, because it only changes at each interaction
    pos_vx = p%v%x > 0._dp ! whether photon is moving in the +ve x direction
    pos_vy = p%v%y > 0._dp ! whether photon is moving in the +ve y direction
    pos_vz = p%v%z > 0._dp ! whether photon is moving in the +ve z direction

    ! Store inv_v to go faster
    if(pos_vx) then
       tx = ( fab%w1(p%icell%i1+1) - p%r%x ) / p%v%x
    else if(p%v%x < 0._dp) then
       tx = ( fab%w1(p%icell%i1)   - p%r%x ) / p%v%x
    else
       tx = huge(1._dp)
    end if

    if(pos_vy) then
       ty = ( fab%w2(p%icell%i2+1) - p%r%y ) / p%v%y
    else if(p%v%y < 0._dp) then
       ty = ( fab%w2(p%icell%i2)   - p%r%y ) / p%v%y
    else
       ty = huge(1._dp)
    end if

    if(pos_vz) then
       tz = ( fab%w3(p%icell%i3+1) - p%r%z ) / p%v%z
    else if(p%v%z < 0._dp) then
       tz = ( fab%w3(p%icell%i3)   - p%r%z ) / p%v%z
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
