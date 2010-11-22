! This approximation is not valid for cells which have emission in them, but
! that's fine, because if they had emission in them we wouldn't be calling this
! here

module grid_pda

  use core_lib
  use type_grid_cell
  use grid_physics
  use grid_geometry
  use dust_main
  use grid_io
  use grid_pda_geometry

  implicit none
  save

  private
  public :: solve_pda

  real(dp), parameter :: tolerance_iter = 1.e-4 ! temperature calculation convergence criterion
  real(dp), parameter :: tolerance_exact = 1.e-5 ! temperature calculation convergence criterion
  real(dp), parameter :: threshold_pda = 0.005 ! maximum number of photons required to use PDA
contains

  subroutine solve_pda()

    implicit none

    real(dp), allocatable :: temperature_prev(:,:)
    logical,allocatable :: do_pda(:)
    real(dp) :: maxdiff
    real(dp) :: mean_n_photons
    real(dp) :: tolerance

    integer :: ic
    integer :: ipda
    type(grid_cell), allocatable :: pda_cells(:)
    integer,allocatable :: id_pda_cell(:)

    mean_n_photons = sum(n_photons) / size(n_photons)

    allocate(temperature_prev(geo%n_cells, n_dust))
    allocate(do_pda(geo%n_cells))

    do_pda = n_photons < max(30,ceiling(threshold_pda*mean_n_photons)) .and. sum(density, dim=2) > 0.

    call check_allowed_pda(do_pda)

    if(.not.any(do_pda)) then
       write(*,'(" [pda] not necessary for this iteration")')
       return
    end if

    if(count(do_pda) < 10000) then
       write(*,'(" [pda] less than 10,000 PDA cells - using Gauss pivot method")')
       tolerance = tolerance_exact
    else
       write(*,'(" [pda] more than 10,000 PDA cells - using iterative method")')
       tolerance = tolerance_iter
    end if

    ! Precompute the cells where the PDA will be computed
    allocate(pda_cells(count(do_pda)))

    allocate(id_pda_cell(geo%n_cells))
    id_pda_cell = -1

    ipda = 0
    do ic=1,geo%n_cells
       if(do_pda(ic)) then
          ipda = ipda + 1
          pda_cells(ipda) = new_grid_cell(ic, geo)
          id_pda_cell(ic) = ipda
       end if
    end do

    do

       temperature_prev = temperature

       call update_mean_temperature()
       call update_dtau_rosseland()

       if(count(do_pda) < 10000) then
          call solve_pda_indiv_exact(pda_cells, id_pda_cell)
       else
          call solve_pda_indiv_iterative(pda_cells)
       end if

       maxdiff = maxval(abs(temperature - temperature_prev) / temperature_prev)

       write(*,'(" [pda] maximum temperature difference: ", ES9.2)') maxdiff

       if(maxdiff < tolerance) exit

    end do

    write(*,'(" [pda] converged")')

    deallocate(do_pda)
    deallocate(temperature_prev)
    
    call update_energy_abs() ! Update energy absorbed in each cell to match temperature

  end subroutine solve_pda

  subroutine solve_pda_indiv_exact(pda_cells, id_pda_cell)

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    integer,intent(in) :: id_pda_cell(:)
    ! Which cells should be used for the PDA

    integer :: direction, wall
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next

    integer :: id
    real(dp) :: coefficient
    real(dp),allocatable :: a(:,:), b(:)
    integer :: id_curr, id_next

    allocate(a(size(pda_cells), size(pda_cells)), b(size(pda_cells)))
    a = 0._dp
    b = 0._dp

    do id_curr=1,size(pda_cells)

       curr = pda_cells(id_curr)

       do wall = 1, geo%n_dim * 2

          direction = int((wall+1)/2)

          next = next_cell(curr, wall)

          dtau_ross_curr = dtau_rosseland(curr%ic, direction)
          dtau_ross_next = dtau_rosseland(next%ic, direction)

          coefficient = 1. / (dtau_ross_curr + dtau_ross_next) / geo%width(curr%ic, direction)
          coefficient = coefficient * geometrical_factor(wall, curr)

          a(id_curr, id_curr) = a(id_curr, id_curr) - coefficient

          if(id_pda_cell(next%ic) > 0) then
             id_next = id_pda_cell(next%ic)
             a(id_next, id_curr) = coefficient
          else
             b(id_curr) = b(id_curr) - coefficient * temperature_mean(next%ic) ** 4.
          end if

       end do

    end do

    call lineq_gausselim(a, b)

    do id_curr=1,size(pda_cells)
       do id=1,n_dust
          if(density(pda_cells(id_curr)%ic, id) > 0._dp) then
             temperature(pda_cells(id_curr)%ic, id) = sqrt(sqrt(b(id_curr)))
          end if
       end do
    end do

    deallocate(a, b)

  end subroutine solve_pda_indiv_exact

  subroutine solve_pda_indiv_iterative(pda_cells)

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    ! Which cells should be used for the PDA

    integer :: direction, wall
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next

    real(dp) :: coefficient
    real(dp) :: a, b
    integer :: id, id_curr

    real(dp),allocatable :: t4mean(:)
    real(dp) :: maxt4diff, t4diff, t4new

    allocate(t4mean(geo%n_cells))
    t4mean = temperature_mean ** 4

    do

       maxt4diff = 0.
       do id_curr=1,size(pda_cells)

          curr = pda_cells(id_curr)

          a = 0._dp
          b = 0._dp

          do wall = 1, geo%n_dim * 2

             direction = int((wall+1)/2)

             next = next_cell(curr, wall)

             dtau_ross_curr = dtau_rosseland(curr%ic, direction)
             dtau_ross_next = dtau_rosseland(next%ic, direction)

             coefficient = 1. / (dtau_ross_curr + dtau_ross_next) / geo%width(curr%ic, direction)
             coefficient = coefficient * geometrical_factor(wall, curr)

             a = a - coefficient
             b = b - coefficient * t4mean(next%ic)

          end do

          t4new = b/a

          t4diff = abs(t4new - t4mean(curr%ic)) / t4mean(curr%ic)
          if(t4diff > maxt4diff) maxt4diff = t4diff

          t4mean(curr%ic) = t4new

       end do

       if(maxt4diff < tolerance_iter) exit

    end do

    do id=1,n_dust
       where(density(:,id) > 0.)
          temperature(:,id) = sqrt(sqrt(t4mean))
       end where
    end do

  end subroutine solve_pda_indiv_iterative

end module grid_pda
