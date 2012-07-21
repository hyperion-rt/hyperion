module grid_geometry

  use core_lib
  use mpi_core
  use mpi_hdf5_io
  use type_grid_cell
  use grid_io, only : read_grid_3d
  use grid_geometry_specific

  implicit none
  save

  private

  ! Imported from grid-specific module
  public :: grid_geometry_debug
  public :: find_cell
  public :: next_cell
  public :: place_in_cell
  public :: in_correct_cell
  public :: random_position_cell
  public :: find_wall
  public :: distance_to_closest_wall
  public :: setup_grid_geometry
  public :: geo
  public :: escaped
  public :: cell_width

  public :: random_cell
  public :: random_masked_cell

  public :: grid_load_pdf_map
  public :: grid_sample_pdf_map
  public :: grid_sample_pdf_map2

  public :: opposite_wall

contains

  type(wall_id) function opposite_wall(id)
    type(wall_id),intent(in) :: id
    opposite_wall%w1 = -id%w1
    opposite_wall%w2 = -id%w2
    opposite_wall%w3 = -id%w3
  end function opposite_wall

  subroutine grid_load_pdf_map(group, path, pdf)

    implicit none

    integer(hid_t),intent(in) :: group
    character(len=*), intent(in) :: path
    type(pdf_discrete_dp), intent(out) :: pdf
    real(dp), allocatable :: map(:)

    ! Read in map from file
    allocate(map(geo%n_cells))
    call read_grid_3d(group, path, map, geo)

    ! Set up PDF to sample map
    call set_pdf(pdf, map)

  end subroutine grid_load_pdf_map

  subroutine grid_sample_pdf_map(pdf, icell)
    implicit none
    type(pdf_discrete_dp), intent(in) :: pdf
    type(grid_cell),intent(out) :: icell
    icell = new_grid_cell(sample_pdf(pdf), geo)
  end subroutine grid_sample_pdf_map

  subroutine grid_sample_pdf_map2(pdf, icell, prob)

    implicit none

    type(pdf_discrete_dp), intent(in) :: pdf
    type(grid_cell),intent(out) :: icell
    integer :: ic
    real(dp) :: prob
    real(dp) :: xi

    do
       call random(xi)
       ic = ceiling(xi * pdf%n)
       prob = pdf%pdf(ic)
       if(prob > 1.e-100_dp) exit
    end do

    icell = new_grid_cell(ic, geo)

    prob = prob * pdf%n

  end subroutine grid_sample_pdf_map2

  type(grid_cell) function random_cell()
    implicit none
    real(dp) :: xi
    integer :: ic
    call random(xi)
    ic = ceiling(xi * geo%n_cells)
    random_cell = new_grid_cell(ic, geo)
  end function random_cell

  type(grid_cell) function random_masked_cell()
    implicit none
    real(dp) :: xi
    integer :: ic
    call random(xi)
    if (geo%masked) then
       ic = geo%mask_map(ceiling(xi * geo%n_masked))
    else
       ic = ceiling(xi * geo%n_cells)
    end if
    random_masked_cell = new_grid_cell(ic, geo)
  end function random_masked_cell

end module grid_geometry
