module type_grid

  use core_lib

  implicit none
  save

  private

  public :: grid_desc
  type grid_desc
     integer :: n_cells, n_dim, n1, n2, n3
     real(dp), allocatable :: w1(:), w2(:), w3(:)
     integer :: start_id
     real(dp) :: volume, area(6), width(3)
     integer, allocatable :: goto_grid(:,:,:)
     integer, allocatable :: goto_level(:,:,:)
     real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  end type grid_desc

  public :: level_desc
  type level_desc
     type(grid_desc),allocatable :: grids(:)
  end type level_desc

  public :: grid_geometry_desc
  type grid_geometry_desc

     character(len=32) :: id
     integer :: n_cells, n_dim
     real(dp), allocatable :: volume(:)
     type(level_desc),allocatable :: levels(:)
     character(len=10) :: type
     real(dp) :: eps

     ! Masking
     logical :: masked = .true.
     integer :: n_masked
     logical, allocatable :: mask(:)
     integer, allocatable :: mask_map(:)

  end type grid_geometry_desc

end module type_grid
