module grid_escape

  ! Estimate, for every cell, the minimum column density (for each dust type)
  ! that a photon has to cross to escape the grid.
  !
  ! Rays are started from random positions in randomly chosen cells, with
  ! random isotropic directions, and are integrated all the way to the edge of
  ! the grid, recording each cell crossed and the column density accumulated
  ! up to the wall through which the ray leaves it. The column density from
  ! that exit wall to the edge of the grid (the total column density along the
  ! ray minus the column density accumulated up to the wall) is a sample of
  ! the column density to escape from that cell. Since column density can only
  ! accumulate along a ray, the minimum over a cell is always found on the
  ! cell walls, so this is where the samples are taken, and every ray provides
  ! a sample for every cell it crosses. Each sample is the exact column density
  ! along a real path, so the minimum found can only be larger than or equal
  ! to the true minimum over the (closed) cell, and converges to it from above
  ! as the number of rays increases. The minimum is kept separately for each
  ! dust type, so that multiplying by the extinction coefficient of each dust
  ! type and summing gives an estimate of the minimum optical depth along any
  ! ray leaving the cell. Rays that intersect a source are discarded, and
  ! cells that do not receive any samples are given a minimum of zero.
  !
  ! With MPI, rays are handed out to processes dynamically using the same
  ! scheduler as for photons, and the minima are then combined so that all
  ! processes have the final values.

  use core_lib, only : dp, idp, error, random
  use mpi_core, only : main_process
  use mpi_routines, only : mp_reset_first, mp_n_photons, mp_allreduce_min_2d
  use type_photon, only : photon
  use type_angle3d, only : random_sphere_angle3d
  use type_vector3d, only : angle3d_to_vector3d, operator(.dot.), operator(+), operator(*)
  use type_grid_cell, only : new_grid_cell, wall_id, no_wall, operator(.eq.)
  use grid_geometry, only : geo, random_position_cell, escaped, find_wall, next_cell, opposite_wall
  use grid_physics, only : density
  use dust_main, only : n_dust
  use sources, only : find_nearest_source

  implicit none
  save

  private
  public :: compute_escape_column_density

  ! Maximum number of cells a ray can cross
  integer, parameter :: max_steps = 1000000

contains

  subroutine compute_escape_column_density(rays_per_cell, column_min)

    implicit none

    ! Average number of rays to start in each cell
    integer, intent(in) :: rays_per_cell

    ! Minimum column density to escape for each cell and dust type
    real(dp), intent(out) :: column_min(:,:)

    integer :: ic, k, n_steps, n_start_cells
    integer(idp) :: n_rays_requested, n_rays_current, n_chunk, ip
    integer, allocatable :: start_cells(:), path_cell(:)
    real(dp), allocatable :: path_column(:,:), total(:)
    type(photon) :: p
    logical :: killed
    real(dp) :: xi

    if(main_process()) write(*,'(" [grid_escape] computing minimum column density to escape each cell")')

    ! Make a list of the cells that rays can start in
    allocate(start_cells(geo%n_cells))
    n_start_cells = 0
    do ic = 1, geo%n_cells
       if(geo%masked) then
          if(.not.geo%mask(ic)) cycle
       end if
       if(geo%volume(ic) == 0._dp) cycle
       n_start_cells = n_start_cells + 1
       start_cells(n_start_cells) = ic
    end do

    allocate(path_cell(max_steps), path_column(n_dust, max_steps), total(n_dust))

    column_min = huge(1._dp)

    n_rays_requested = int(rays_per_cell, idp) * int(n_start_cells, idp)
    n_rays_current = 0
    call mp_reset_first()

    do

       call mp_n_photons(n_rays_requested, n_rays_current, max(n_rays_requested / 100_idp, 1_idp), n_chunk)

       if(n_chunk == 0) exit

       do ip = 1, n_chunk

          call random(xi)
          ic = start_cells(max(ceiling(xi * real(n_start_cells, dp)), 1))

          p%icell = new_grid_cell(ic, geo)
          p%in_cell = .true.
          p%on_wall = .false.
          call random_position_cell(p%icell, p%r)
          call random_sphere_angle3d(p%a)
          call angle3d_to_vector3d(p%a, p%v)

          call trace(p, n_steps, path_cell, path_column, total, killed)

          if(killed) cycle

          ! Every cell crossed (including the starting cell) gets the column
          ! density from the wall through which the ray leaves it
          do k = 1, n_steps
             column_min(path_cell(k), :) = min(column_min(path_cell(k), :), &
                  &                            max(total - path_column(:, k), 0._dp))
          end do

       end do

    end do

    call mp_allreduce_min_2d(column_min)

    ! Cells that did not receive any samples (for example if all the rays
    ! started in them intersected a source) are given a minimum column density
    ! of zero, so that they are never treated as hidden
    where(column_min == huge(1._dp))
       column_min = 0._dp
    end where

  end subroutine compute_escape_column_density

  subroutine trace(p_orig, n_steps, path_cell, path_column, total, killed)

    ! Integrate the column density along a ray to the edge of the grid,
    ! recording each cell crossed and the column density accumulated up to the
    ! wall through which the ray leaves it

    implicit none

    type(photon), intent(in) :: p_orig
    integer, intent(out) :: n_steps
    integer, intent(out) :: path_cell(:)
    real(dp), intent(out) :: path_column(:,:), total(:)
    logical, intent(out) :: killed

    type(photon) :: p
    type(wall_id) :: id_min
    real(dp) :: tmin, t_source
    integer :: source_id
    logical :: radial

    p = p_orig
    killed = .false.
    n_steps = 0
    total = 0._dp

    if(escaped(p)) return

    call find_nearest_source(p%r, p%v, t_source, source_id)
    if(t_source < huge(1._dp)) then
       killed = .true.
       return
    end if

    radial = (p%r .dot. p%v) > 0.

    do

       call find_wall(p, radial, tmin, id_min)

       if(id_min == no_wall) then
          killed = .true.
          return
       end if

       p%r = p%r + tmin * p%v
       total = total + density(p%icell%ic, :) * tmin

       n_steps = n_steps + 1
       if(n_steps > size(path_cell)) call error('trace', 'too many cells crossed')
       path_cell(n_steps) = p%icell%ic
       path_column(:, n_steps) = total

       p%on_wall = .true.
       p%icell = next_cell(p%icell, id_min, intersection=p%r)
       p%on_wall_id = opposite_wall(id_min)

       if(escaped(p)) exit

    end do

  end subroutine trace

end module grid_escape
