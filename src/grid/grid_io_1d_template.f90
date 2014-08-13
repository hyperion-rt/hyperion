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
    grid_exists = mp_path_exists(group, name)
  end function grid_exists

  !!@FOR real(sp):sp real(dp):dp integer:int integer(idp):int8

  subroutine read_grid_4d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array(group, path, array)

    if(any(is_nan(array))) call error("read_grid_4d", "NaN values in 4D array")

  end subroutine read_grid_4d_<T>

  subroutine read_grid_3d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(out) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array(group, path, array)

    if(any(is_nan(array))) call error("read_grid_3d", "NaN values in 3D array")

  end subroutine read_grid_3d_<T>

  subroutine write_grid_4d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, array)
    call mp_write_keyword(group, path, 'geometry', geo%id)

  end subroutine write_grid_4d_<T>

  subroutine write_grid_3d_<T>(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    @T, intent(in) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, array)
    call mp_write_keyword(group,path, 'geometry', geo%id)

  end subroutine write_grid_3d_<T>

  !!@END FOR

end module grid_io
