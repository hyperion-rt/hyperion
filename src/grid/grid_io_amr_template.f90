module grid_io

  use core_lib
  use mpi_hdf5_io
  use type_grid

  implicit none
  save

  private
  public :: grid_exists
  public :: read_grid_3d
  public :: read_grid_4d
  public :: write_grid_3d
  public :: write_grid_4d

  interface read_grid_3d
     module procedure read_grid_3d_sp
     module procedure read_grid_3d_dp
     module procedure read_grid_3d_int
     module procedure read_grid_3d_int8
  end interface read_grid_3d

  interface read_grid_4d
     module procedure read_grid_4d_sp
     module procedure read_grid_4d_dp
     module procedure read_grid_4d_int
     module procedure read_grid_4d_int8
  end interface read_grid_4d

  interface write_grid_3d
     module procedure write_grid_3d_sp
     module procedure write_grid_3d_dp
     module procedure write_grid_3d_int
     module procedure write_grid_3d_int8
  end interface write_grid_3d

  interface write_grid_4d
     module procedure write_grid_4d_sp
     module procedure write_grid_4d_dp
     module procedure write_grid_4d_int
     module procedure write_grid_4d_int8
  end interface write_grid_4d

contains

  logical function grid_exists(group, name)
    implicit none
    integer(hid_t),intent(in) :: group
    character(len=*),intent(in) :: name
    integer(hid_t) :: g_level, g_grid
    if(mp_path_exists(group, 'level_00001')) then
       g_level = mp_open_group(group, 'level_00001')
       if(mp_path_exists(g_level, 'grid_00001')) then
          g_grid = mp_open_group(g_level, 'grid_00001')
          if(mp_path_exists(g_grid, name)) then
             grid_exists = .true.
          else
             grid_exists = .false.
          end if
       else
          grid_exists = .false.
       end if
    else
       grid_exists = .false.
    end if
  end function grid_exists

  !!@FOR real(sp):sp real(dp):dp integer:int integer(idp):int8

  subroutine read_grid_4d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    @T, allocatable :: array4d(:,:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, igrid, idust
    type(level_desc), pointer :: level
    type(grid_desc), pointer :: grid

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          write(full_path, '("level_", I5.5, "/grid_", I5.5,"/")') ilevel, igrid
          full_path = trim(full_path)//trim(path)
          call mp_read_array_auto(group, full_path, array4d)
          if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")
          do idust=1,size(array4d, 4)
             where(grid%goto_grid(1:grid%n1,1:grid%n2,1:grid%n3) > 0)
                array4d(:,:,:,idust) = 0
             end where
          end do
          array(grid%start_id:grid%start_id + grid%n_cells - 1, :) = reshape(array4d, (/grid%n_cells, size(array, 2)/))
       end do
    end do

    ! The following three lines provide a workaround for the PGI Fortran
    ! compiler, which otherwise crashes with the following error:
    ! 0: RESHAPE: result type != SOURCE type
  contains
    subroutine test()
    end subroutine test

  end subroutine read_grid_4d_<T>

  subroutine read_grid_3d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(out) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    @T, allocatable :: array3d(:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, igrid
    type(level_desc), pointer :: level
    type(grid_desc), pointer :: grid

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          write(full_path, '("level_", I5.5, "/grid_ ", I5.5,"/")') ilevel, igrid
          full_path = trim(full_path)//trim(path)
          call mp_read_array_auto(group, full_path, array3d)
          if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")
          where(grid%goto_grid(1:grid%n1,1:grid%n2,1:grid%n3) > 0)
             array3d(:,:,:) = 0
          end where
          array(grid%start_id:grid%start_id + grid%n_cells - 1) = reshape(array3d, (/grid%n_cells/))
       end do
    end do

    ! The following three lines provide a workaround for the PGI Fortran
    ! compiler, which otherwise crashes with the following error:
    ! 0: RESHAPE: result type != SOURCE type
  contains
    subroutine test()
    end subroutine test

  end subroutine read_grid_3d_<T>

  subroutine write_grid_4d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_grid
    character(len=100) :: name
    integer :: ilevel, igrid
    type(level_desc), pointer :: level
    type(grid_desc), pointer :: grid

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("level_", I5.5)') ilevel
       if(mp_path_exists(group, name)) then
          g_level = mp_open_group(group, name)
       else
          g_level = mp_create_group(group, name)
       end if
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          write(name, '("grid_", I5.5)') igrid
          if(mp_path_exists(g_level, name)) then
             g_grid = mp_open_group(g_level, name)
          else
             g_grid = mp_create_group(g_level, name)
          end if
          call mp_write_array(g_grid, path, reshape(array(grid%start_id:grid%start_id + grid%n_cells - 1, :), &
               &                                     (/grid%n1, grid%n2, grid%n3, size(array,2)/)))
          call mp_close_group(g_grid)
       end do
       call mp_close_group(g_level)

    end do

  end subroutine write_grid_4d_<T>

  subroutine write_grid_3d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(in) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_grid
    character(len=100) :: name
    integer :: ilevel, igrid
    type(level_desc), pointer :: level
    type(grid_desc), pointer :: grid

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("level_", I5.5)') ilevel
       if(mp_path_exists(group, name)) then
          g_level = mp_open_group(group, name)
       else
          g_level = mp_create_group(group, name)
       end if
       do igrid=1,size(level%grids)
          grid => level%grids(igrid)
          write(name, '("grid_", I5.5)') igrid
          if(mp_path_exists(g_level, name)) then
             g_grid = mp_open_group(g_level, name)
          else
             g_grid = mp_create_group(g_level, name)
          end if
          call mp_write_array(g_grid, path, reshape(array(grid%start_id:grid%start_id + grid%n_cells - 1), (/grid%n1, grid%n2, grid%n3/)))
          call mp_close_group(g_grid)
       end do
       call mp_close_group(g_level)

    end do

  end subroutine write_grid_3d_<T>

  !!@END FOR

end module grid_io
