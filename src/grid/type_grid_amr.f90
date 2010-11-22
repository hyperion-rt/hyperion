module type_grid

  use core_lib

  implicit none
  save

  private

  public :: fab_desc
  type fab_desc
     integer :: n_cells, n_dim, n1, n2, n3
     real(dp), allocatable :: w1(:), w2(:), w3(:)
     integer :: start_id
     real(dp) :: volume, area(6), width(3)
     integer, allocatable :: goto_fab(:,:,:)
     integer, allocatable :: goto_level(:,:,:)
     real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  end type fab_desc

  public :: level_desc
  type level_desc
     type(fab_desc),allocatable :: fabs(:)
  end type level_desc

  public :: grid_geometry_desc
  type grid_geometry_desc
     character(len=32) :: id
     integer :: n_cells, n_dim
     real(dp), allocatable :: volume(:)
     real(dp), allocatable :: area(:, :)
     real(dp), allocatable :: width(:, :)
     type(level_desc),allocatable :: levels(:)
     character(len=10) :: type
     real(dp) :: eps
  end type grid_geometry_desc

end module type_grid
