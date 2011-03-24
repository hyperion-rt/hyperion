module mpi_io

  use core_lib
  use mpi
  use mpi_core

  implicit none
  save

  private

  integer :: ierr

  public :: mp_read_keyword
  interface mp_read_keyword
     module procedure mp_read_keyword_mpi_real4
     module procedure mp_read_keyword_mpi_real8
     module procedure mp_read_keyword_mpi_integer4
     module procedure mp_read_keyword_mpi_integer8
     module procedure mp_read_keyword_mpi_logical
     module procedure mp_read_keyword_mpi_character
  end interface mp_read_keyword

  public :: mp_table_read_column_auto
  interface mp_table_read_column_auto
     module procedure mp_table_read_column_auto_1d_mpi_real4
     module procedure mp_table_read_column_auto_1d_mpi_real8
     module procedure mp_table_read_column_auto_1d_mpi_integer4
     module procedure mp_table_read_column_auto_1d_mpi_integer8
     module procedure mp_table_read_column_auto_1d_mpi_character
     module procedure mp_table_read_column_auto_2d_mpi_real4
     module procedure mp_table_read_column_auto_2d_mpi_real8
     module procedure mp_table_read_column_auto_2d_mpi_integer4
     module procedure mp_table_read_column_auto_2d_mpi_integer8
  end interface mp_table_read_column_auto

  public :: mp_read_array_auto
  interface mp_read_array_auto
     module procedure mp_read_array_auto_2d_mpi_real4
     module procedure mp_read_array_auto_2d_mpi_real8
     module procedure mp_read_array_auto_2d_mpi_integer4
     module procedure mp_read_array_auto_2d_mpi_integer8
     module procedure mp_read_array_auto_3d_mpi_real4
     module procedure mp_read_array_auto_3d_mpi_real8
     module procedure mp_read_array_auto_3d_mpi_integer4
     module procedure mp_read_array_auto_3d_mpi_integer8
     module procedure mp_read_array_auto_4d_mpi_real4
     module procedure mp_read_array_auto_4d_mpi_real8
     module procedure mp_read_array_auto_4d_mpi_integer4
     module procedure mp_read_array_auto_4d_mpi_integer8
     module procedure mp_read_array_auto_5d_mpi_real4
     module procedure mp_read_array_auto_5d_mpi_real8
     module procedure mp_read_array_auto_5d_mpi_integer4
     module procedure mp_read_array_auto_5d_mpi_integer8
     module procedure mp_read_array_auto_6d_mpi_real4
     module procedure mp_read_array_auto_6d_mpi_real8
     module procedure mp_read_array_auto_6d_mpi_integer4
     module procedure mp_read_array_auto_6d_mpi_integer8
  end interface mp_read_array_auto

contains

  subroutine mp_read_keyword_mpi_character(handle, path, name, value)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    character(len=*),intent(out) :: value
    if(main_process()) call hdf5_read_keyword(handle, path, name, value)
    call mpi_bcast(value, len(value), mpi_character, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_keyword_mpi_character

  subroutine mp_read_keyword_mpi_logical(handle, path, name, value)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    logical,intent(out) :: value
    if(main_process()) call hdf5_read_keyword(handle, path, name, value)
    call mpi_bcast(value, 1, mpi_logical, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_keyword_mpi_logical

  subroutine mp_table_read_column_auto_1d_mpi_character(handle, path, name, array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    character(len=*),allocatable,intent(out) :: array(:)
    integer :: n1
    if(main_process()) call hdf5_table_read_column_auto(handle, path, name, array)
    n1 = size(array)
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1))
    call mpi_bcast(array, n1, mpi_character, rank_main, mpi_comm_world, ierr)
  end subroutine mp_table_read_column_auto_1d_mpi_character

  !!@FOR real(sp):mpi_real4 real(dp):mpi_real8 integer:mpi_integer4 integer(idp):mpi_integer8

  subroutine mp_read_keyword_<T>(handle, path, name, value)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    @T,intent(out) :: value
    if(main_process()) call hdf5_read_keyword(handle, path, name, value)
    call mpi_bcast(value, 1, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_keyword_<T>

  subroutine mp_table_read_column_auto_1d_<T>(handle, path, name, array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    @T,allocatable,intent(out) :: array(:)
    integer :: n1
    if(main_process()) call hdf5_table_read_column_auto(handle, path, name, array)
    n1 = size(array)
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1))
    call mpi_bcast(array, n1, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_table_read_column_auto_1d_<T>

  subroutine mp_table_read_column_auto_2d_<T>(handle, path, name, array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path, name
    @T,allocatable,intent(out) :: array(:, :)
    integer :: n1, n2
    if(main_process()) then
       call hdf5_table_read_column_auto(handle, path, name, array)
       n1 = size(array, 1)
       n2 = size(array, 2)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2))
    call mpi_bcast(array, n1*n2, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_table_read_column_auto_2d_<T>

  subroutine mp_read_array_auto_2d_<T>(handle,path,array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path
    @T,allocatable,intent(out) :: array(:, :)
    integer :: n1, n2
    if(main_process()) then
       call hdf5_read_array_auto(handle,path,array)
       n1 = size(array, 1)
       n2 = size(array, 2)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2))
    call mpi_bcast(array, n1*n2, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_array_auto_2d_<T>

  subroutine mp_read_array_auto_3d_<T>(handle,path,array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path
    @T,allocatable,intent(out) :: array(:, :, :)
    integer :: n1, n2, n3
    if(main_process()) then
       call hdf5_read_array_auto(handle,path,array)
       n1 = size(array, 1)
       n2 = size(array, 2)
       n3 = size(array, 3)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n3, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2, n3))
    call mpi_bcast(array, n1*n2*n3, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_array_auto_3d_<T>

  subroutine mp_read_array_auto_4d_<T>(handle,path,array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path
    @T,allocatable,intent(out) :: array(:, :, :, :)
    integer :: n1, n2, n3, n4
    if(main_process()) then
       call hdf5_read_array_auto(handle,path,array)
       n1 = size(array, 1)
       n2 = size(array, 2)
       n3 = size(array, 3)
       n4 = size(array, 4)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n3, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n4, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2, n3, n4))
    call mpi_bcast(array, n1*n2*n3*n4, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_array_auto_4d_<T>

  subroutine mp_read_array_auto_5d_<T>(handle,path,array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path
    @T,allocatable,intent(out) :: array(:, :, :, :, :)
    integer :: n1, n2, n3, n4, n5
    if(main_process()) then
       call hdf5_read_array_auto(handle,path,array)
       n1 = size(array, 1)
       n2 = size(array, 2)
       n3 = size(array, 3)
       n4 = size(array, 4)
       n5 = size(array, 4)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n3, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n4, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n5, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2, n3, n4, n5))
    call mpi_bcast(array, n1*n2*n3*n4*n5, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_array_auto_5d_<T>

  subroutine mp_read_array_auto_6d_<T>(handle,path,array)
    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path
    @T,allocatable,intent(out) :: array(:, :, :, :, :, :)
    integer :: n1, n2, n3, n4, n5, n6
    if(main_process()) then
       call hdf5_read_array_auto(handle,path,array)
       n1 = size(array, 1)
       n2 = size(array, 2)
       n3 = size(array, 3)
       n4 = size(array, 4)
       n5 = size(array, 5)
       n6 = size(array, 6)
    end if
    call mpi_bcast(n1, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n2, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n3, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n4, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n5, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(n6, 1, mpi_integer4, rank_main, mpi_comm_world, ierr)
    if(.not. main_process()) allocate(array(n1, n2, n3, n4, n5, n6))
    call mpi_bcast(array, n1*n2*n3*n4*n5*n6, <T>, rank_main, mpi_comm_world, ierr)
  end subroutine mp_read_array_auto_6d_<T>


  !!@END FOR

end module mpi_io

