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

  public :: setup_grid_geometry
  type(grid_geometry_desc),public,target :: geo

contains

  ! AMR Helper functions

  real(dp) function cell_width(cell, idir)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: idir
    cell_width = geo%levels(cell%ilevel)%grids(cell%igrid)%width(idir)
  end function cell_width

  logical function grids_intersect(grid1, grid2)
    ! Given two grids, do they intersect anywhere?
    implicit none
    type(grid_desc),intent(in) :: grid1, grid2
    grids_intersect = .false.
    if(grid1%xmax < grid2%xmin) return
    if(grid1%xmin > grid2%xmax) return
    if(grid1%ymax < grid2%ymin) return
    if(grid1%ymin > grid2%ymax) return
    if(grid1%zmax < grid2%zmin) return
    if(grid1%zmin > grid2%zmax) return
    grids_intersect = .true.
  end function grids_intersect

  logical function grids_close(grid1, grid2)
    ! Given two grids, are they within one half cell from each other
    type(grid_desc),intent(in) :: grid1, grid2
    grids_close = .false.
    if(grid1%xmax < grid2%xmin - grid2%width(1)*0.5_dp) return
    if(grid1%xmin > grid2%xmax + grid2%width(1)*0.5_dp) return
    if(grid1%ymax < grid2%ymin - grid2%width(2)*0.5_dp) return
    if(grid1%ymin > grid2%ymax + grid2%width(2)*0.5_dp) return
    if(grid1%zmax < grid2%zmin - grid2%width(3)*0.5_dp) return
    if(grid1%zmin > grid2%zmax + grid2%width(3)*0.5_dp) return
    grids_close = .true.
  end function grids_close

  logical function in_grid(grid, r)
    ! Given a grid and a position r, determines if the position lies inside the grid.
    ! Could be optimized by finding distance from center (abs(x)-x0) < dx, although not clear which is faster.
    implicit none
    type(grid_desc),intent(in) :: grid
    type(vector3d_dp),intent(in) :: r
    in_grid = .false.
    if(r%x < grid%xmin) return
    if(r%x > grid%xmax) return
    if(r%y < grid%ymin) return
    if(r%y > grid%ymax) return
    if(r%z < grid%zmin) return
    if(r%z > grid%zmax) return
    in_grid = .true.
  end function in_grid

  integer function locate_grid(level, r) result(igrid)
    ! Given a level and a position r, determines if the position lies in
    ! any of the grids, and if so returns the grid ID in the level. If the
    ! position does not lie in any grid, the function returns -1.
    implicit none
    type(level_desc),intent(in) :: level
    type(vector3d_dp),intent(in) :: r
    do igrid=1,size(level%grids)
       if(in_grid(level%grids(igrid), r)) return
    end do
    igrid = -1
  end function locate_grid

  subroutine read_grid(group, grid)

    ! Given and HDF5 group containing a grid, this subroutine reads and returns it

    implicit none

    integer(hid_t),intent(in) :: group
    type(grid_desc),intent(out) :: grid

    call mp_read_keyword(group, '.', 'n1', grid%n1)
    call mp_read_keyword(group, '.', 'n2', grid%n2)
    call mp_read_keyword(group, '.', 'n3', grid%n3)
    grid%n_cells = grid%n1 * grid%n2 * grid%n3
    grid%n_dim = 3

    allocate(grid%w1(grid%n1+1))
    call mp_read_keyword(group, '.', 'xmin', grid%xmin)
    call mp_read_keyword(group, '.', 'xmax', grid%xmax)
    call linspace(grid%xmin, grid%xmax, grid%w1)

    allocate(grid%w2(grid%n2+1))
    call mp_read_keyword(group, '.', 'ymin', grid%ymin)
    call mp_read_keyword(group, '.', 'ymax', grid%ymax)
    call linspace(grid%ymin, grid%ymax, grid%w2)

    allocate(grid%w3(grid%n3+1))
    call mp_read_keyword(group, '.', 'zmin', grid%zmin)
    call mp_read_keyword(group, '.', 'zmax', grid%zmax)
    call linspace(grid%zmin, grid%zmax, grid%w3)

    grid%width(1) = (grid%xmax-grid%xmin)/real(grid%n1, dp)
    grid%width(2) = (grid%ymax-grid%ymin)/real(grid%n2, dp)
    grid%width(3) = (grid%zmax-grid%zmin)/real(grid%n3, dp)

    grid%area(1:2) = grid%width(2) * grid%width(3)
    grid%area(3:4) = grid%width(1) * grid%width(3)
    grid%area(5:6) = grid%width(1) * grid%width(2)

    grid%volume = grid%width(1) * grid%width(2) * grid%width(3)

    allocate(grid%goto_grid(0:grid%n1+1, 0:grid%n2+1, 0:grid%n3+1))
    allocate(grid%goto_level(0:grid%n1+1, 0:grid%n2+1, 0:grid%n3+1))

    grid%goto_grid = 0
    grid%goto_level = 0

  end subroutine read_grid

  subroutine read_level(group, level)

    ! Given and HDF5 group containing a level, this subroutine reads and returns it

    integer(hid_t),intent(in) :: group
    type(level_desc),intent(out) :: level
    integer :: ngrids, igrid
    integer(hid_t) :: g_grid
    character(len=100) :: grid_name

    call mp_read_keyword(group, '.', 'ngrids', ngrids)

    allocate(level%grids(ngrids))

    do igrid=1,ngrids
       write(grid_name, '("grid_", I5.5)') igrid
       g_grid = mp_open_group(group, grid_name)
       call read_grid(g_grid, level%grids(igrid))
       call mp_close_group(g_grid)
    end do

  end subroutine read_level

  ! Standard Grid Methods

  logical function aligned(x1, x2, dx)
    implicit none
    real(dp),intent(in) :: x1, x2, dx
    real(dp) :: r
    r = mod(abs(x1 - x2), dx)
    if(r > 0.5_dp * dx) r = dx - r
    aligned = abs(r / dx) < 1.e-8
  end function aligned

  subroutine setup_grid_geometry(group)

    implicit none

    integer(hid_t),intent(in) :: group

    integer :: nlevels, ilevel, ilevel1, ilevel2
    integer(hid_t) :: g_level
    character(len=100) :: level_name

    integer :: igrid, igrid1, igrid2, icell, iv

    integer :: start_id
    type(level_desc), pointer :: level, level1, level2
    type(grid_desc), pointer :: grid, grid1, grid2, grid_ref
    integer :: i1, i2, i3
    type(vector3d_dp) :: r
    real(dp) :: min_width
    real(dp) :: ref(3)

    type(grid_cell) :: cell

    character(len=1000) :: message

    ! Read geometry file
    call mp_read_keyword(group, '.', "geometry", geo%id)
    call mp_read_keyword(group, '.', "grid_type", geo%type)

    if(trim(geo%type).ne.'amr') call error("setup_grid_geometry","grid is not AMR cartesian")

    if(main_process()) write(*,'(" [setup_grid_geometry] Reading AMR cartesian grid")')

    ! Read the number of levels
    call mp_read_keyword(group, '.', 'nlevels', nlevels)

    ! Allocate the levels
    allocate(geo%levels(nlevels))

    ! Loop through the levels and read all the grids
    do ilevel=1,nlevels
       write(level_name, '("level_", I5.5)') ilevel
       g_level = mp_open_group(group, level_name)
       call read_level(g_level, geo%levels(ilevel))
       call mp_close_group(g_level)
    end do

    ! Check that all grids in a given level have a common width and have edges
    ! that line up on a grid
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       grid_ref => level%grids(1)
       do igrid=2,size(level%grids)
          grid => level%grids(igrid)
          if(abs(grid%width(1) - grid_ref%width(1)) > 1.e-10 * grid%width(1)) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have differing cell widths in the x direction (", ES11.4, " and ", ES11.4, " respectively)")') 1, igrid, ilevel, grid_ref%width(1), grid%width(1)
             call error("setup_grid_geometry", trim(message))
          end if
          if(abs(grid%width(2) - grid_ref%width(2)) > 1.e-10 * grid%width(2)) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have differing cell widths in the y direction (", ES11.4, " and ", ES11.4, " respectively)")') 1, igrid, ilevel, grid_ref%width(2), grid%width(2)
             call error("setup_grid_geometry", trim(message))
          end if
          if(abs(grid%width(3) - grid_ref%width(3)) > 1.e-10 * grid%width(3)) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have differing cell widths in the z direction (", ES11.4, " and ", ES11.4, " respectively)")') 1, igrid, ilevel, grid_ref%width(3), grid%width(3)
             call error("setup_grid_geometry", trim(message))
          end if
          if(.not.aligned(grid%xmin, grid_ref%xmin, grid_ref%width(1))) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have edges that are not separated by an integer number of cells in the x direction")') 1, igrid, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
          if(.not.aligned(grid%ymin, grid_ref%ymin, grid_ref%width(2))) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have edges that are not separated by an integer number of cells in the y direction")') 1, igrid, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
          if(.not.aligned(grid%zmin, grid_ref%zmin, grid_ref%width(3))) then
             write(message, '("Grids ", I0, " and ", I0, " in level ", I0, " have edges that are not separated by an integer number of cells in the z direction")') 1, igrid, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
       end do
    end do

    ! Check that the refinement factor in each direction is an integer
    do ilevel=1,size(geo%levels) - 1
       level1 => geo%levels(ilevel)
       level2 => geo%levels(ilevel + 1)
       grid1 => level1%grids(1)
       grid2 => level2%grids(1)
       ref = grid1%width / grid2%width
       if(abs(ref(1) - nint(ref(1))) > 1.e-10) then
          write(message, '("Refinement factor in the x direction between level ", I0, " and level ", I0, " is not an integer (",F0.3,")")') ilevel, ilevel + 1, ref(1)
          call error("setup_grid_geometry", trim(message))
       end if
       if(abs(ref(2) - nint(ref(2))) > 1.e-10) then
          write(message, '("Refinement factor in the y direction between level ", I0, " and level ", I0, " is not an integer (",F0.3,")")') ilevel, ilevel + 1, ref(2)
          call error("setup_grid_geometry", trim(message))
       end if
       if(abs(ref(3) - nint(ref(3))) > 1.e-10) then
          write(message, '("Refinement factor in the z direction between level ", I0, " and level ", I0, " is not an integer (",F0.3,")")') ilevel, ilevel + 1, ref(3)
          call error("setup_grid_geometry", trim(message))
       end if
    end do

    ! Check that grid edges line up with cells in parent level
    do ilevel=1,size(geo%levels) - 1
       level1 => geo%levels(ilevel)
       level2 => geo%levels(ilevel + 1)
       grid_ref => level1%grids(1)
       do igrid=1,size(level2%grids)
          grid => level2%grids(igrid)
          if(.not.aligned(grid%xmin, grid_ref%xmin, grid_ref%width(1))) then
             write(message, '("Grid ", I0, " in level ", I0, " is not aligned with cells in level ", I0, " in the x direction")') igrid, ilevel + 1, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
          if(.not.aligned(grid%ymin, grid_ref%ymin, grid_ref%width(2))) then
             write(message, '("Grid ", I0, " in level ", I0, " is not aligned with cells in level ", I0, " in the y direction")') igrid, ilevel + 1, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
          if(.not.aligned(grid%zmin, grid_ref%zmin, grid_ref%width(3))) then
             write(message, '("Grid ", I0, " in level ", I0, " is not aligned with cells in level ", I0, " in the z direction")') igrid, ilevel + 1, ilevel
             call error("setup_grid_geometry", trim(message))
          end if
       end do
    end do

    ! Find the unique ID of the first cell in each grid
    start_id = 1
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          grid%start_id = start_id
          start_id = start_id + grid%n_cells
       end do
    end do

    ! Set the total number of cells
    geo%n_cells = start_id - 1

    ! Compute volume for all cells
    allocate(geo%volume(geo%n_cells))
    min_width = huge(1._dp)
    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          if(any(grid%width < min_width)) min_width = minval(grid%width)
          do icell=grid%start_id,grid%start_id+grid%n_cells-1
             geo%volume(icell) = grid%volume
             if(geo%volume(icell)==0.) stop
          end do
       end do
    end do
    if(any(geo%volume==0._dp)) call error('setup_grid_geometry','all volumes should be greater than zero')

    ! Set number of dimensions (is this needed by PDA?)
    geo%n_dim = 3

    ! Set the step size for stepping out of grids. Choose half the smallest width in the grid.
    geo%eps = min_width / 2._dp

    ! Find the level ID, grid ID, and internal 3D coordinates for each unique ID
    call preset_cell_id(geo)

    ! For each grid, find all cells in level below that overlap, and set flag in level below to indicate this

    do ilevel1=size(geo%levels)-1, 1, -1
       level1 => geo%levels(ilevel1)
       level2 => geo%levels(ilevel1+1)
       do igrid1=1,size(level1%grids)
          grid1 => level1%grids(igrid1)
          do igrid2=1,size(level2%grids)
             grid2 => level2%grids(igrid2)
             if(grids_intersect(grid1, grid2)) then
                do i1=1,grid1%n1
                   r%x = 0.5_dp * (grid1%w1(i1) + grid1%w1(i1+1))
                   do i2=1,grid1%n2
                      r%y = 0.5_dp * (grid1%w2(i2) + grid1%w2(i2+1))
                      do i3=1,grid1%n3
                         r%z = 0.5_dp * (grid1%w3(i3) + grid1%w3(i3+1))
                         if(in_grid(grid2, r)) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel1+1
                         end if
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

    ! Now for all neighboring cells, find out what grid they end up in if they go one step outside the grid

    do ilevel1 = 1, size(geo%levels)
       level1 => geo%levels(ilevel1)
       do igrid1=1,size(level1%grids)
          grid1 => level1%grids(igrid1)

          do ilevel2 = ilevel1, 1, -1
             level2 => geo%levels(ilevel2)
             do igrid2=1,size(level2%grids)
                grid2 => level2%grids(igrid2)

                if(grids_close(grid1, grid2).and.(igrid1.ne.igrid2 .or. ilevel1.ne.ilevel2)) then

                   ! xmin side
                   i1 = 0
                   r%x = grid1%xmin - grid1%width(1) * 0.5_dp
                   do i2=1,grid1%n2
                      r%y = 0.5_dp * (grid1%w2(i2) + grid1%w2(i2+1))
                      do i3=1,grid1%n3
                         r%z = 0.5_dp * (grid1%w3(i3) + grid1%w3(i3+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! xmax side
                   i1 = grid1%n1+1
                   r%x = grid1%xmax + grid1%width(1) * 0.5_dp
                   do i2=1,grid1%n2
                      r%y = 0.5_dp * (grid1%w2(i2) + grid1%w2(i2+1))
                      do i3=1,grid1%n3
                         r%z = 0.5_dp * (grid1%w3(i3) + grid1%w3(i3+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymin side
                   i2 = 0
                   r%y = grid1%ymin - grid1%width(2) * 0.5_dp
                   do i1=1,grid1%n1
                      r%x = 0.5_dp * (grid1%w1(i1) + grid1%w1(i1+1))
                      do i3=1,grid1%n3
                         r%z = 0.5_dp * (grid1%w3(i3) + grid1%w3(i3+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymax side
                   i2 = grid1%n2+1
                   r%y = grid1%ymax + grid1%width(2) * 0.5_dp
                   do i1=1,grid1%n1
                      r%x = 0.5_dp * (grid1%w1(i1) + grid1%w1(i1+1))
                      do i3=1,grid1%n3
                         r%z = 0.5_dp * (grid1%w3(i3) + grid1%w3(i3+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! zmin side
                   i3 = 0
                   r%z = grid1%zmin - grid1%width(3) * 0.5_dp
                   do i1=1,grid1%n1
                      r%x = 0.5_dp * (grid1%w1(i1) + grid1%w1(i1+1))
                      do i2=1,grid1%n2
                         r%y = 0.5_dp * (grid1%w2(i2) + grid1%w2(i2+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                   ! ymax side
                   i3 = grid1%n3+1
                   r%z = grid1%zmax + grid1%width(3) * 0.5_dp
                   do i1=1,grid1%n1
                      r%x = 0.5_dp * (grid1%w1(i1) + grid1%w1(i1+1))
                      do i2=1,grid1%n2
                         r%y = 0.5_dp * (grid1%w2(i2) + grid1%w2(i2+1))
                         if(in_grid(grid2, r) .and. grid1%goto_grid(i1, i2, i3) == 0) then
                            grid1%goto_grid(i1, i2, i3) = igrid2
                            grid1%goto_level(i1, i2, i3) = ilevel2
                         end if
                      end do
                   end do

                end if

             end do
          end do

       end do
    end do

    ! Figure out which cells are valid, and construct an index to map valid
    ! cell IDs to original IDs
    geo%masked = .true.
    allocate(geo%mask(geo%n_cells))
    do icell=1,geo%n_cells
       cell = new_grid_cell(icell, geo)
       grid => geo%levels(cell%ilevel)%grids(cell%igrid)
       geo%mask(icell) = grid%goto_grid(cell%i1, cell%i2, cell%i3) == 0
    end do
    geo%n_masked = count(geo%mask)
    allocate(geo%mask_map(geo%n_masked))
    iv = 0
    do icell=1,geo%n_cells
       if(geo%mask(icell)) then
          iv = iv + 1
          geo%mask_map(iv) = icell
       end if
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

  recursive type(grid_cell) function find_position_in_grid(r, ilevel, igrid) result(cell)
    implicit none
    type(vector3d_dp),intent(in) :: r
    integer,intent(in) :: ilevel, igrid
    integer :: i1, i2, i3, ilevel_new, igrid_new
    type(grid_desc), pointer :: grid
    grid => geo%levels(ilevel)%grids(igrid)
    i1 = ipos2(grid%xmin, grid%xmax, r%x, grid%n1)
    i2 = ipos2(grid%ymin, grid%ymax, r%y, grid%n2)
    i3 = ipos2(grid%zmin, grid%zmax, r%z, grid%n3)
    ilevel_new = grid%goto_level(i1, i2, i3)
    igrid_new = grid%goto_grid(i1, i2, i3)
    if(ilevel_new == 0) then
       cell = new_grid_cell(i1, i2, i3, ilevel, igrid, geo)
    else
       cell = find_position_in_grid(r, ilevel_new, igrid_new)
    end if
  end function find_position_in_grid

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
    integer :: igrid
    if(debug) write(*,'(" [debug] find_cell")')
    igrid = locate_grid(geo%levels(1), r) ! Find grid in level 1
    icell = find_position_in_grid(r, 1, igrid) ! Refine position
  end function find_cell_position

  subroutine place_in_cell(p)
    implicit none
    type(photon),intent(inout) :: p
    p%icell = find_cell(p)
    if(p%icell == invalid_cell.or.p%icell == outside_cell) then
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
    escaped_cell = cell == outside_cell
  end function escaped_cell

  type(grid_cell) function next_cell_int(cell, direction, intersection)
    implicit none
    type(grid_cell),intent(in) :: cell
    integer,intent(in) :: direction
    type(vector3d_dp),optional,intent(in) :: intersection
    type(vector3d_dp) :: r

    integer :: i1, i2, i3, ilevel, igrid
    type(grid_desc), pointer :: grid

    grid => geo%levels(cell%ilevel)%grids(cell%igrid)

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

    if(grid%goto_level(i1, i2, i3) == 0) then
       if(i1==0.or.i1==grid%n1+1.or.i2==0.or.i2==grid%n2+1.or.i3==0.or.i3==grid%n3+1) then
          next_cell_int = outside_cell
       else
          next_cell_int = new_grid_cell(i1, i2, i3, cell%ilevel, cell%igrid, geo)
       end if
    else
       ilevel = grid%goto_level(i1, i2, i3)
       igrid = grid%goto_grid(i1, i2, i3)
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
       next_cell_int = find_position_in_grid(r, ilevel, igrid)
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
    ! Numerical issue here with photons that are currently on walls. Don't try and find right level and icell, assume they are correct.

    implicit none
    type(photon),intent(in) :: p
    type(grid_cell) :: icell_actual
    real(dp) :: frac
    type(grid_desc), pointer :: grid
    grid => geo%levels(p%icell%ilevel)%grids(p%icell%igrid)
    icell_actual = p%icell
    icell_actual%i1 = ipos(grid%xmin, grid%xmax, p%r%x, grid%n1)
    icell_actual%i2 = ipos(grid%ymin, grid%ymax, p%r%y, grid%n2)
    icell_actual%i3 = ipos(grid%zmin, grid%zmax, p%r%z, grid%n3)
    if(p%on_wall) then

       in_correct_cell = .true.

       if(p%on_wall_id%w1 == -1) then
          frac = (p%r%x - grid%w1(p%icell%i1)) / (grid%w1(p%icell%i1+1) - grid%w1(p%icell%i1))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
       else if(p%on_wall_id%w1 == +1) then
          frac = (p%r%x - grid%w1(p%icell%i1+1)) / (grid%w1(p%icell%i1+1) - grid%w1(p%icell%i1))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i1 == p%icell%i1
       end if

       if(p%on_wall_id%w2 == -1) then
          frac = (p%r%y - grid%w2(p%icell%i2)) / (grid%w2(p%icell%i2+1) - grid%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
       else if(p%on_wall_id%w2 == +1) then
          frac = (p%r%y - grid%w2(p%icell%i2+1)) / (grid%w2(p%icell%i2+1) - grid%w2(p%icell%i2))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
       else
          in_correct_cell = in_correct_cell .and. icell_actual%i2 == p%icell%i2
       end if

       if(p%on_wall_id%w3 == -1) then
          frac = (p%r%z - grid%w3(p%icell%i3)) / (grid%w3(p%icell%i3+1) - grid%w3(p%icell%i3))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
       else if(p%on_wall_id%w3 == +1) then
          frac = (p%r%z - grid%w3(p%icell%i3+1)) / (grid%w3(p%icell%i3+1) - grid%w3(p%icell%i3))
          in_correct_cell = in_correct_cell .and. abs(frac) < 1.e-3_dp
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
    type(grid_desc), pointer :: grid
    grid => geo%levels(icell%ilevel)%grids(icell%igrid)
    call random(x)
    call random(y)
    call random(z)
    pos%x = x * (grid%w1(icell%i1+1) - grid%w1(icell%i1)) + grid%w1(icell%i1)
    pos%y = y * (grid%w2(icell%i2+1) - grid%w2(icell%i2)) + grid%w2(icell%i2)
    pos%z = z * (grid%w3(icell%i3+1) - grid%w3(icell%i3)) + grid%w3(icell%i3)
  end subroutine random_position_cell

  real(dp) function distance_to_closest_wall(p) result(d)

    implicit none

    type(photon),intent(in) :: p

    real(dp) :: d1,d2,d3,d4,d5,d6

    type(grid_desc), pointer :: grid

    grid => geo%levels(p%icell%ilevel)%grids(p%icell%igrid)

    d1 = p%r%x - grid%w1(p%icell%i1)
    d2 = grid%w1(p%icell%i1+1) - p%r%x
    d3 = p%r%y - grid%w2(p%icell%i2)
    d4 = grid%w2(p%icell%i2+1) - p%r%y
    d5 = p%r%z - grid%w3(p%icell%i3)
    d6 = grid%w3(p%icell%i3+1) - p%r%z

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

    type(grid_desc), pointer :: grid

    grid => geo%levels(p%icell%ilevel)%grids(p%icell%igrid)

    ! Can store this in photon, because it only changes at each interaction
    pos_vx = p%v%x > 0._dp ! whether photon is moving in the +ve x direction
    pos_vy = p%v%y > 0._dp ! whether photon is moving in the +ve y direction
    pos_vz = p%v%z > 0._dp ! whether photon is moving in the +ve z direction

    ! Store inv_v to go faster
    if(pos_vx) then
       tx = ( grid%w1(p%icell%i1+1) - p%r%x ) / p%v%x
    else if(p%v%x < 0._dp) then
       tx = ( grid%w1(p%icell%i1)   - p%r%x ) / p%v%x
    else
       tx = huge(1._dp)
    end if

    if(pos_vy) then
       ty = ( grid%w2(p%icell%i2+1) - p%r%y ) / p%v%y
    else if(p%v%y < 0._dp) then
       ty = ( grid%w2(p%icell%i2)   - p%r%y ) / p%v%y
    else
       ty = huge(1._dp)
    end if

    if(pos_vz) then
       tz = ( grid%w3(p%icell%i3+1) - p%r%z ) / p%v%z
    else if(p%v%z < 0._dp) then
       tz = ( grid%w3(p%icell%i3)   - p%r%z ) / p%v%z
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

  end subroutine find_wall

end module grid_geometry_specific
