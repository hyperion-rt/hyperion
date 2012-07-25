module type_grid

  use core_lib

  implicit none
  save

  private

  public :: grid_geometry_desc
  type grid_geometry_desc
     character(len=32) :: id
     integer :: n_cells, n_dim, n1, n2, n3
     real(dp), allocatable :: volume(:)
     real(dp), allocatable :: w1(:), w2(:), w3(:)
     real(dp), allocatable :: ew1(:), ew2(:), ew3(:)
     real(dp), allocatable :: wr2(:), wcost(:), wsint(:), wtant(:), wtant2(:), wtanp(:)
     real(dp), allocatable :: r(:), dr(:), dr2(:), dr3(:), t(:), dt(:), dcost(:), dphi(:)
     integer :: midplane = -1
     character(len=10) :: type

     ! Masking
     logical :: masked = .false.
     integer :: n_masked
     logical, allocatable :: mask(:)
     integer, allocatable :: mask_map(:)

  end type grid_geometry_desc

end module type_grid
