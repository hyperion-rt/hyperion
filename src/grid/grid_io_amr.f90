! MD5 of template: 5eb8f79665a237db77deb18f58d27ade
module grid_io

  use core_lib
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
    integer(hid_t) :: g_level, g_fab
    if(hdf5_path_exists(group, 'Level 1')) then
       g_level = hdf5_open_group(group, 'Level 1')
       if(hdf5_path_exists(g_level, 'Fab 1')) then
          g_fab = hdf5_open_group(g_level, 'Fab 1')
          if(hdf5_path_exists(g_fab, name)) then
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


  subroutine read_grid_4d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(idp), allocatable :: array4d(:,:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array4d)
          if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1, :) = reshape(array4d, (/fab%n_cells, size(array, 2)/))
       end do
    end do

  end subroutine read_grid_4d_int8

  subroutine read_grid_3d_int8(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(idp), allocatable :: array3d(:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array3d)
          if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1) = reshape(array3d, (/fab%n_cells/))
       end do
    end do

  end subroutine read_grid_3d_int8

  subroutine write_grid_4d_int8(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1, :), (/fab%n1, fab%n2, fab%n3, size(array,2)/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_4d_int8

  subroutine write_grid_3d_int8(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer(idp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1), (/fab%n1, fab%n2, fab%n3/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_3d_int8


  subroutine read_grid_4d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer, allocatable :: array4d(:,:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array4d)
          if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1, :) = reshape(array4d, (/fab%n_cells, size(array, 2)/))
       end do
    end do

  end subroutine read_grid_4d_int

  subroutine read_grid_3d_int(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(out) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer, allocatable :: array3d(:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array3d)
          if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1) = reshape(array3d, (/fab%n_cells/))
       end do
    end do

  end subroutine read_grid_3d_int

  subroutine write_grid_4d_int(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1, :), (/fab%n1, fab%n2, fab%n3, size(array,2)/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_4d_int

  subroutine write_grid_3d_int(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    integer, intent(in) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1), (/fab%n1, fab%n2, fab%n3/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_3d_int


  subroutine read_grid_4d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    real(dp), allocatable :: array4d(:,:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array4d)
          if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1, :) = reshape(array4d, (/fab%n_cells, size(array, 2)/))
       end do
    end do

  end subroutine read_grid_4d_dp

  subroutine read_grid_3d_dp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    real(dp), allocatable :: array3d(:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array3d)
          if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1) = reshape(array3d, (/fab%n_cells/))
       end do
    end do

  end subroutine read_grid_3d_dp

  subroutine write_grid_4d_dp(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1, :), (/fab%n1, fab%n2, fab%n3, size(array,2)/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_4d_dp

  subroutine write_grid_3d_dp(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1), (/fab%n1, fab%n2, fab%n3/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_3d_dp


  subroutine read_grid_4d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(out) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    real(sp), allocatable :: array4d(:,:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array4d)
          if(any(is_nan(array4d))) call error("read_grid_4d", "NaN values in 4D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1, :) = reshape(array4d, (/fab%n_cells, size(array, 2)/))
       end do
    end do

  end subroutine read_grid_4d_sp

  subroutine read_grid_3d_sp(group, path, array, geo)

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(out) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    real(sp), allocatable :: array3d(:,:,:)
    character(len=100) :: full_path
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(full_path, '("Level ", I0, "/Fab ", I0,"/")') ilevel, ifab
          full_path = trim(full_path)//trim(path)
          call hdf5_read_array_auto(group, full_path, array3d)
          if(any(is_nan(array3d))) call error("read_grid_3d", "NaN values in 3D array")
          array(fab%start_id:fab%start_id + fab%n_cells - 1) = reshape(array3d, (/fab%n_cells/))
       end do
    end do

  end subroutine read_grid_3d_sp

  subroutine write_grid_4d_sp(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(in) :: array(:,:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1, :), (/fab%n1, fab%n2, fab%n3, size(array,2)/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_4d_sp

  subroutine write_grid_3d_sp(group, path, array, geo) 

    implicit none

    integer(hid_t), intent(in) :: group
    character(len=*), intent(in) :: path
    real(sp), intent(in) :: array(:)
    type(grid_geometry_desc),intent(in),target :: geo
    integer(hid_t) :: g_level, g_fab
    character(len=100) :: name
    integer :: ilevel, ifab
    type(level_desc), pointer :: level
    type(fab_desc), pointer :: fab

    do ilevel=1,size(geo%levels)
       level => geo%levels(ilevel)
       write(name, '("Level ", I0)') ilevel
       if(hdf5_path_exists(group, name)) then
          g_level = hdf5_open_group(group, name)
       else
          g_level = hdf5_create_group(group, name)
       end if
       do ifab=1,size(level%fabs)
          fab => level%fabs(ifab)
          write(name, '("Fab ", I0)') ifab
          if(hdf5_path_exists(g_level, name)) then
             g_fab = hdf5_open_group(g_level, name)
          else
             g_fab = hdf5_create_group(g_level, name)
          end if
          call hdf5_write_array(g_fab, path, reshape(array(fab%start_id:fab%start_id + fab%n_cells - 1), (/fab%n1, fab%n2, fab%n3/)))
          call hdf5_close_group(g_fab)
       end do
       call hdf5_close_group(g_level)

    end do

  end subroutine write_grid_3d_sp


end module grid_io
