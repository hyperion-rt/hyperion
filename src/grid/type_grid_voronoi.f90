module type_grid

  use core_lib

  implicit none
  save

  private

  ! Define refinable cell
  public :: cell
  type cell
     real(dp) :: x,y,z  ! position of the point defining the cell
     integer,allocatable,dimension(:) :: neighbors
  end type cell

  ! Define Octree descriptor
  public :: grid_geometry_desc
  type grid_geometry_desc

     character(len=32) :: id
     integer :: n_cells, n_dim
     type(cell),allocatable,dimension(:) :: cells
     real(dp), allocatable :: volume(:)
     character(len=10) :: type
     real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
     
  end type grid_geometry_desc

end module type_grid

