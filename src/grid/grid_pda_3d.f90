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

  real(dp), parameter :: tolerance_iter = 1.e-4 ! energy calculation convergence criterion
  real(dp), parameter :: tolerance_exact = 1.e-5 ! energy calculation convergence criterion
  real(dp), parameter :: threshold_pda = 0.005 ! maximum number of photons required to use PDA

  real(dp), allocatable :: e_mean(:)

contains

  real(dp) elemental function difference_ratio(a, b)
    implicit none
    real(dp), intent(in) :: a, b
    difference_ratio = max(a/b, b/a)
  end function difference_ratio

  subroutine update_specific_energy(ic)
    implicit none
    integer,intent(in) :: ic
    integer :: id
    real(dp) :: s_prev, s, smin, smax
    do id=1,n_dust

       s = specific_energy(ic, id)

       smin = d(id)%specific_energy(1)
       smax = d(id)%specific_energy(d(id)%n_e)

       if(e_mean(ic) < smin / kappa_planck(id, smin)) then
          call warn("update_specific_energy", "specific energy in PDA below minimum allowed by dust type - resetting")
          s = smin
       else if (e_mean(ic) > smax / kappa_planck(id, smax)) then
          call warn("update_specific_energy", "specific energy in PDA above maximum allowed by dust type - resetting")
          s = smax
       else
          do
             s_prev = s
             s = e_mean(ic) * kappa_planck(id, s)
             if(difference_ratio(s, s_prev) - 1._dp < 1.e-5_dp) exit
          end do
       end if
       specific_energy(ic, id) = s
    end do
  end subroutine update_specific_energy

  subroutine update_e_mean(ic)
    implicit none
    integer,intent(in) :: ic
    integer :: id
    e_mean(ic) = 0.
    if(sum(density(ic, :)) > 0._dp) then
       do id=1,n_dust
          e_mean(ic) = e_mean(ic) + density(ic,id) * specific_energy(ic,id) / kappa_planck(id, specific_energy(ic,id))
       end do
       e_mean(ic) = e_mean(ic) / sum(density(ic, :))
    end if
  end subroutine update_e_mean

  subroutine solve_pda()

    implicit none

    real(dp),allocatable :: specific_energy_prev(:,:)
    logical,allocatable :: do_pda(:)
    real(dp) :: maxdiff
    real(dp) :: mean_n_photons
    real(dp) :: tolerance

    integer :: ic
    integer :: ipda
    type(grid_cell), allocatable :: pda_cells(:)
    integer,allocatable :: id_pda_cell(:)

    mean_n_photons = sum(n_photons) / size(n_photons)

    allocate(specific_energy_prev(geo%n_cells, n_dust))
    allocate(do_pda(geo%n_cells))

    do_pda = n_photons < max(30,ceiling(threshold_pda*mean_n_photons)) .and. sum(density, dim=2) > 0._dp

    call check_allowed_pda(do_pda)

    if(.not.any(do_pda)) then
       write(*,'(" [pda] not necessary for this iteration")')
       return
    end if

    if(count(do_pda) < 10000) then
       write(*,'(" [pda] fewer than 10,000 PDA cells - using Gauss pivot method")')
       tolerance = tolerance_exact
    else
       write(*,'(" [pda] more than 10,000 PDA cells - using iterative method")')
       tolerance = tolerance_iter
    end if

    allocate(e_mean(geo%n_cells))

    do ic=1,geo%n_cells
       call update_e_mean(ic)
    end do

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

    specific_energy_prev = specific_energy

    do

       specific_energy_prev = specific_energy

       if(count(do_pda) < 10000) then
          call solve_pda_indiv_exact(pda_cells, id_pda_cell)
       else
          call solve_pda_indiv_iterative(pda_cells)
       end if

       maxdiff = maxval(abs(specific_energy - specific_energy_prev) / specific_energy_prev)

       write(*,'(" [pda] maximum energy difference: ", ES9.2)') maxdiff

       if(maxdiff < tolerance) exit

    end do

    write(*,'(" [pda] converged")')

    deallocate(do_pda)
    deallocate(specific_energy_prev)
    deallocate(e_mean)

    call update_energy_abs_tot()

    call check_energy_abs()

  end subroutine solve_pda

  real(dp) function dtau_rosseland(cell, idir)
    implicit none
    type(grid_cell), intent(in) :: cell
    integer,intent(in) :: idir
    integer :: id
    dtau_rosseland = 0._dp
    do id=1,n_dust
       dtau_rosseland = dtau_rosseland + density(cell%ic,id) * chi_rosseland(id, specific_energy(cell%ic,id)) * cell_width(cell,idir)
    end do
  end function dtau_rosseland

  subroutine solve_pda_indiv_exact(pda_cells, id_pda_cell)

    implicit none

    type(grid_cell),intent(in) :: pda_cells(:)
    integer,intent(in) :: id_pda_cell(:)
    ! Which cells should be used for the PDA

    integer :: direction, wall
    type(grid_cell) :: curr, next
    real(dp) :: dtau_ross_curr, dtau_ross_next, dtau_sum

    integer :: ic
    real(dp) :: coefficient
    real(dp),allocatable :: a(:,:), b(:)
    integer :: id_curr, id_next

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_e_mean(ic)
    end do

    allocate(a(size(pda_cells), size(pda_cells)), b(size(pda_cells)))
    a = 0._dp
    b = 0._dp

    do id_curr=1,size(pda_cells)

       curr = pda_cells(id_curr)

       do wall = 1, geo%n_dim * 2

          direction = int((wall+1)/2)

          next = next_cell(curr, wall)

          dtau_ross_curr = dtau_rosseland(curr, direction)
          dtau_ross_next = dtau_rosseland(next, direction)

          dtau_sum = dtau_ross_curr + dtau_ross_next

          ! If the optical depth is too small, we have to reset it to avoid
          ! issues.
          if(dtau_sum < 1e-100_dp) dtau_sum = 1e-100_dp

          coefficient = 1. / dtau_sum / cell_width(curr, direction)
          coefficient = coefficient * geometrical_factor(wall, curr)

          a(id_curr, id_curr) = a(id_curr, id_curr) - coefficient

          if(id_pda_cell(next%ic) > 0) then
             id_next = id_pda_cell(next%ic)
             a(id_next, id_curr) = coefficient
          else
             b(id_curr) = b(id_curr) - coefficient * e_mean(next%ic)
          end if

       end do

    end do

    call lineq_gausselim(a, b)

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       e_mean(ic) = b(id_curr)
       call update_specific_energy(ic)
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
    integer :: id_curr, ic

    real(dp) :: max_e_diff, e_diff, e_new

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_e_mean(ic)
    end do

    do

       max_e_diff = 0.
       do id_curr=1,size(pda_cells)

          curr = pda_cells(id_curr)

          a = 0._dp
          b = 0._dp

          do wall = 1, geo%n_dim * 2

             direction = int((wall+1)/2)

             next = next_cell(curr, wall)

             dtau_ross_curr = dtau_rosseland(curr, direction)
             dtau_ross_next = dtau_rosseland(next, direction)

             coefficient = 1. / (dtau_ross_curr + dtau_ross_next) / cell_width(curr, direction)
             coefficient = coefficient * geometrical_factor(wall, curr)

             a = a - coefficient
             b = b - coefficient * e_mean(next%ic)

          end do

          e_new = b/a

          e_diff = abs(e_new - e_mean(curr%ic)) / e_mean(curr%ic)
          if(e_diff > max_e_diff) max_e_diff = e_diff

          e_mean(curr%ic) = e_new

       end do

       if(max_e_diff < tolerance_iter) exit

    end do

    do id_curr=1,size(pda_cells)
       ic = pda_cells(id_curr)%ic
       call update_specific_energy(ic)
    end do

  end subroutine solve_pda_indiv_iterative

end module grid_pda
