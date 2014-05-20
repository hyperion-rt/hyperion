module type_grid

  use kdtree2_module
  use core_lib

  implicit none
  save

  private

  ! Define refinable cell
  public :: cell
  type cell
      integer :: id
      type(vector3d_dp) :: r
     integer,allocatable,dimension(:) :: neighbors

     ! Bounding box
     real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax

  end type cell

  ! Define Octree descriptor
  public :: grid_geometry_desc
  type grid_geometry_desc

     type(kdtree2), pointer :: tree
     character(len=32) :: id
     integer :: n_cells, n_dim
     type(cell),allocatable,dimension(:) :: cells
     real(dp), allocatable :: volume(:)
     character(len=10) :: type
     real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
     
     ! Masking
     logical :: masked = .false.
     integer :: n_masked
     logical, allocatable :: mask(:)
     integer, allocatable :: mask_map(:)

  end type grid_geometry_desc

end module type_grid

