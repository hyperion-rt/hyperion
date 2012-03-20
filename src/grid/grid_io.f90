! MD5 of template: aaf3fd4ae68c1556a0ea47b490321b92
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


  subroutine read_grid_4d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo
    integer(idp), allocatable :: array4d(:,:,:,:)
    integer :: n_cells, n_dust

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array4d)

    if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")

    n_cells = size(array, 1)
    n_dust = size(array, 2)

    array = reshape(array4d, (/n_cells, n_dust/))

  end subroutine read_grid_4d_int8

  subroutine read_grid_3d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo
    integer(idp), allocatable :: array3d(:, :, :)
    integer :: n_cells

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array3d)

    if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")

    n_cells = size(array)

    array = reshape(array3d, (/n_cells/))

  end subroutine read_grid_3d_int8

  subroutine write_grid_4d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, reshape(array, (/geo%n1, geo%n2, geo%n3, size(array,2)/)))
    call mp_write_keyword(group, path, 'geometry', geo%id)

  end subroutine write_grid_4d_int8

  subroutine write_grid_3d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group,path, reshape(array, (/geo%n1, geo%n2, geo%n3/)))
    call mp_write_keyword(group,path, 'geometry', geo%id)

  end subroutine write_grid_3d_int8


  subroutine read_grid_4d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo
    integer, allocatable :: array4d(:,:,:,:)
    integer :: n_cells, n_dust

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array4d)

    if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")

    n_cells = size(array, 1)
    n_dust = size(array, 2)

    array = reshape(array4d, (/n_cells, n_dust/))

  end subroutine read_grid_4d_int

  subroutine read_grid_3d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(out) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo
    integer, allocatable :: array3d(:, :, :)
    integer :: n_cells

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array3d)

    if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")

    n_cells = size(array)

    array = reshape(array3d, (/n_cells/))

  end subroutine read_grid_3d_int

  subroutine write_grid_4d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, reshape(array, (/geo%n1, geo%n2, geo%n3, size(array,2)/)))
    call mp_write_keyword(group, path, 'geometry', geo%id)

  end subroutine write_grid_4d_int

  subroutine write_grid_3d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(in) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group,path, reshape(array, (/geo%n1, geo%n2, geo%n3/)))
    call mp_write_keyword(group,path, 'geometry', geo%id)

  end subroutine write_grid_3d_int


  subroutine read_grid_4d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo
    real(dp), allocatable :: array4d(:,:,:,:)
    integer :: n_cells, n_dust

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array4d)

    if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")

    n_cells = size(array, 1)
    n_dust = size(array, 2)

    array = reshape(array4d, (/n_cells, n_dust/))

  end subroutine read_grid_4d_dp

  subroutine read_grid_3d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo
    real(dp), allocatable :: array3d(:, :, :)
    integer :: n_cells

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array3d)

    if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")

    n_cells = size(array)

    array = reshape(array3d, (/n_cells/))

  end subroutine read_grid_3d_dp

  subroutine write_grid_4d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, reshape(array, (/geo%n1, geo%n2, geo%n3, size(array,2)/)))
    call mp_write_keyword(group, path, 'geometry', geo%id)

  end subroutine write_grid_4d_dp

  subroutine write_grid_3d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group,path, reshape(array, (/geo%n1, geo%n2, geo%n3/)))
    call mp_write_keyword(group,path, 'geometry', geo%id)

  end subroutine write_grid_3d_dp


  subroutine read_grid_4d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo
    real(sp), allocatable :: array4d(:,:,:,:)
    integer :: n_cells, n_dust

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array4d)

    if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")

    n_cells = size(array, 1)
    n_dust = size(array, 2)

    array = reshape(array4d, (/n_cells, n_dust/))

  end subroutine read_grid_4d_sp

  subroutine read_grid_3d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo
    real(sp), allocatable :: array3d(:, :, :)
    integer :: n_cells

    character(len=32) :: geometry_id_check

    call mp_read_keyword(group,path, 'geometry', geometry_id_check)
    if(geometry_id_check.ne.geo%id) then
       call error("read_grid", "geometry IDs do not match")
    end if
    call mp_read_array_auto(group,path, array3d)

    if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")

    n_cells = size(array)

    array = reshape(array3d, (/n_cells/))

  end subroutine read_grid_3d_sp

  subroutine write_grid_4d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group, path, reshape(array, (/geo%n1, geo%n2, geo%n3, size(array,2)/)))
    call mp_write_keyword(group, path, 'geometry', geo%id)

  end subroutine write_grid_4d_sp

  subroutine write_grid_3d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in) :: geo

    call mp_write_array(group,path, reshape(array, (/geo%n1, geo%n2, geo%n3/)))
    call mp_write_keyword(group,path, 'geometry', geo%id)

  end subroutine write_grid_3d_sp


end module grid_io
