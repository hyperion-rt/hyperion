module mpi_routines

  use mpi
  use mpi_core
  use core_lib
  use grid_physics
  use setup
  use performance

  implicit none
  save

  private

  ! MPI errors
  integer :: ierr

  ! MPI status
  integer,dimension(mpi_status_size) :: status

  ! Variable to signal first-time use for mp_n_photons
  logical :: first = .true.

  ! Pre-defined tags
  integer,parameter :: tag1 = 50, tag2 = 60

  logical :: debug = .false.

  public :: mp_reset_first
  public :: mp_n_photons
  public :: mp_collect_physical_arrays
  public :: mp_broadcast_specific_energy
  public :: mp_collect_images
  public :: mp_broadcast_convergence
  public :: mp_set_random_seed

  real(dp) :: time_curr
  integer(idp) :: n_completed, n_photons_chunk
  integer(idp) :: n_steps = 10
  integer(idp) :: n_stats_last

  public :: mp_sync
  interface mp_sync
     module procedure mp_sync_integer4
     module procedure mp_sync_integer8
     module procedure mp_sync_real4
     module procedure mp_sync_real8
  end interface mp_sync

contains

  subroutine mp_reset_first()
    implicit none
    first = .true.
    time_curr = 0._dp
    n_completed = 0
    n_photons_chunk = 0
    n_stats_last = 0
  end subroutine mp_reset_first

  subroutine mp_n_photons(n_photons_tot, n_photons_curr, n_photons_stats, n_photons)

    implicit none

    ! Total number of photons, and size of a chunk
    integer(idp),intent(in) :: n_photons_tot, n_photons_stats

    ! Number of photons requested so far
    integer(idp),intent(inout) :: n_photons_curr

    ! Number of photons to compute this time
    integer(idp),intent(out) :: n_photons

    ! Number of photons used for MPI transfer
    integer(idp),volatile,allocatable :: n_photons_send(:)

    ! Loop variable and dummy variable
    integer :: ir
    real(dp), volatile :: dum_dp

    ! Flag for MPI test
    logical :: flag

    ! Request variables
    integer,allocatable,save :: request(:)
    integer :: request_dum

    ! Whether a process has been asked to stop
    logical,allocatable :: stopped(:)
    logical,allocatable,save :: started(:)

    real(dp),save :: time1 = -1._dp
    real(dp) :: time2
    real(dp),allocatable,volatile,save :: dtime(:)

    if(debug) write(*,'("[mpi_routines] rank ",I0," requesting number of photons")') rank

    if(time1 < 0._dp) call cpu_time(time1)
    call cpu_time(time2)

    select case(rank)
    case(rank_main)

       ! Find the number of photons per chunk
       if(n_photons_chunk == 0) then
          n_photons_chunk = max(nint(real(n_photons_tot, dp) / real(nproc, dp) / real(n_steps, dp)), 1)
       end if

       time_curr = time_curr + time2-time1

       if(.not.allocated(request)) allocate(request(nproc-1))
       if(.not.allocated(dtime)) allocate(dtime(nproc-1))
       if(.not.allocated(n_photons_send)) allocate(n_photons_send(nproc-1))

       if(.not.allocated(stopped)) then
          allocate(stopped(nproc-1))
          stopped = .false.
       end if

       if(.not.allocated(started)) then
          allocate(started(nproc-1))
          started = .false.
       end if

       ! Loop over processes
       do ir=1,nproc-1

          if(.not.first) then
             ! Test whether previous request was received
             call mpi_test(request(ir), flag, status, ierr)
             if(flag) time_curr = time_curr + dtime(ir)
          else
             ! Set flag to true to force sending a request
             flag = .true.
          end if

          if (flag) then

             ! If all photons have been requested, we exit the loop and wrap things up
             if(n_photons_curr == n_photons_tot) exit
             if(n_photons_curr > n_photons_tot) stop "n_photons_curr > n_photons_tot"

             ! Find how many photons to request and increment counter
             if(n_photons_tot - n_photons_curr <= n_photons_chunk * nproc) then
                n_photons_send(ir) = max(nint(real(n_photons_chunk, dp) / 10._dp), 1)
             else
                n_photons_send(ir) = n_photons_chunk
             end if

             if(n_photons_curr + n_photons_send(ir) > n_photons_tot) n_photons_send(ir) = n_photons_tot - n_photons_curr

             n_photons_curr = n_photons_curr + n_photons_send(ir)

             ! Send number of photons and initialize receive for status check
             call mpi_isend(n_photons_send(ir), 1, mpi_integer8, ir, tag1, mpi_comm_world, request_dum, ierr)
             call mpi_irecv(dtime(ir), 1, mpi_real8, ir, tag2, mpi_comm_world, request(ir), ierr)

             if(first) started(ir) = .true.

          end if

       end do

       if(n_photons_curr > n_photons_tot) stop "n_photons_curr > n_photons_tot"

       if(n_photons_curr == n_photons_tot) then

          if(debug) write(*,'("[mpi_routines] master rank now checking all clients are stopped")')

          ! Loop until all processes are finished
          do while(.not.all(stopped))

             call microsleep(200000)

             ! Loop over processes
             do ir=1,nproc-1

                if(started(ir).and..not.stopped(ir)) then

                   ! Test if last valid request was received
                   call mpi_test(request(ir), flag, status, ierr)

                   ! If it was, send out a final request
                   if(flag) then
                      n_photons_send(ir) = 0_idp
                      call mpi_isend(n_photons_send(ir), 1, mpi_integer8, ir, tag1, mpi_comm_world, request_dum, ierr)
                      call mpi_irecv(dum_dp, 1, mpi_real8, ir, tag2, mpi_comm_world, request(ir), ierr)
                      stopped(ir) = .true.
                      if(debug) write(*,'("[mpi_routines] send abort signal to rank ",I0)') ir
                   else
                      if(debug) write(*,'("[mpi_routines] rank ",I0," is not ready for abort signal")') ir
                   end if

                else if(.not.started(ir)) then
                   n_photons_send(ir) = 0_idp
                   call mpi_isend(n_photons_send(ir), 1, mpi_integer8, ir, tag1, mpi_comm_world, request_dum, ierr)
                   call mpi_irecv(dum_dp, 1, mpi_real8, ir, tag2, mpi_comm_world, request(ir), ierr)
                   started(ir) = .true.
                   stopped(ir) = .true.
                   if(debug) write(*,'("[mpi_routines] send abort signal to rank ",I0)') ir
                end if


             end do
          end do


          if(debug) write(*,'("[mpi_routines] master rank now waiting for all jobs to complete")')

          ! Wait for all acknowledgements of final request
          call mpi_waitall(nproc-1, request, mpi_statuses_ignore, ierr)

          ! Set number of photons to 0 for main process too
          n_photons = 0

       else

          ! Set number of photons to a fraction of the normal chunk size
          n_photons = max(nint(real(n_photons_chunk, dp) / 10._dp), 1)

          if(n_photons_curr + n_photons > n_photons_tot) n_photons = n_photons_tot - n_photons_curr

          n_photons_curr = n_photons_curr + n_photons

       end if

       if(first) first=.false.

       if(n_photons_stats > 0) then
          if(n_photons_curr >= n_stats_last + n_photons_stats) then
             if(n_photons_curr > 0) call perf_numbers(n_photons_curr, time_curr)
             n_stats_last = n_photons_curr - mod(n_photons_curr, n_photons_stats)
          end if
       else
          if(n_photons_curr >= n_stats_last + n_photons_chunk) then
             if(n_photons_curr > 0) call perf_numbers(n_photons_curr, time_curr)
             n_stats_last = n_photons_curr - mod(n_photons_curr, n_photons_chunk)
          end if
       end if

    case default

       ! Receive number of photons and send acknowledgments
       call mpi_recv(n_photons, 1, mpi_integer8, rank_main, tag1, mpi_comm_world, status, ierr)
       call mpi_isend(time2-time1, 1, mpi_real8, rank_main, tag2, mpi_comm_world, request_dum, ierr)

       if(n_photons > n_photons_tot) stop "n_photons > n_photons_tot"

    end select

    time1 = time2

    if(debug) write(*,'("[mpi_routines] rank ",I0," will compute ",I0," photons")') rank,n_photons

  end subroutine mp_n_photons

  subroutine mp_set_random_seed(seed)
    implicit none
    integer :: seed
    call set_seed(seed + rank)
  end subroutine mp_set_random_seed

  subroutine mp_collect_physical_arrays()

    implicit none
    real(dp) :: tmp
    real(dp) :: dummy_dp
    integer(idp) :: dummy_idp
    real(dp),allocatable :: tmp_2d(:,:)
    integer(idp),allocatable :: tmp_int_1d(:)

    if(main_process()) then
       allocate(tmp_2d(size(specific_energy_sum,1),size(specific_energy_sum,2)))
       call mpi_reduce(specific_energy_sum, tmp_2d, size(specific_energy_sum), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
       specific_energy_sum = tmp_2d
       deallocate(tmp_2d)
    else
       call mpi_reduce(specific_energy_sum, dummy_dp, size(specific_energy_sum), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
    end if

    if(allocated(n_photons)) then
       if(main_process()) then
          allocate(tmp_int_1d(size(n_photons,1)))
          call mpi_reduce(n_photons, tmp_int_1d, size(n_photons), mpi_integer8, mpi_sum, rank_main, mpi_comm_world, ierr)
          n_photons = tmp_int_1d
          deallocate(tmp_int_1d)
       else
          call mpi_reduce(n_photons, dummy_idp, size(n_photons), mpi_integer8, mpi_sum, rank_main, mpi_comm_world, ierr)
       end if
    end if

  end subroutine mp_collect_physical_arrays

  subroutine mp_broadcast_specific_energy()

    implicit none

    call mpi_bcast(specific_energy, size(specific_energy), mpi_real8, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(energy_abs_tot, size(energy_abs_tot), mpi_real8, rank_main, mpi_comm_world, ierr)

  end subroutine mp_broadcast_specific_energy

  subroutine mp_broadcast_convergence(converged)
    implicit none
    logical,intent(inout) :: converged
    call mpi_bcast(converged, 1, mpi_logical, rank_main, mpi_comm_world, ierr)
  end subroutine mp_broadcast_convergence

  subroutine mp_sync_integer4(value)
    implicit none
    integer,intent(inout) :: value
    integer :: tmp
    call mpi_allreduce(value, tmp, 1, mpi_integer4, mpi_sum, mpi_comm_world, ierr)
    value = tmp
  end subroutine mp_sync_integer4

  subroutine mp_sync_integer8(value)
    implicit none
    integer(idp),intent(inout) :: value
    integer(idp) :: tmp
    call mpi_allreduce(value, tmp, 1, mpi_integer8, mpi_sum, mpi_comm_world, ierr)
    value = tmp
  end subroutine mp_sync_integer8

  subroutine mp_sync_real4(value)
    implicit none
    real(sp),intent(inout) :: value
    real(sp) :: tmp
    call mpi_allreduce(value, tmp, 1, mpi_real4, mpi_sum, mpi_comm_world, ierr)
    value = tmp
  end subroutine mp_sync_real4

  subroutine mp_sync_real8(value)
    implicit none
    real(dp),intent(inout) :: value
    real(dp) :: tmp
    call mpi_allreduce(value, tmp, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    value = tmp
  end subroutine mp_sync_real8

  subroutine mp_collect_images()

    use binned_images
    use peeled_images

    implicit none

    integer :: ip

    real(dp),allocatable :: cube5d(:,:,:,:,:)
    real(dp),allocatable :: cube6d(:,:,:,:,:,:)

    if(make_binned_images) then

       if(binned_image%compute_sed) then

          allocate(cube5d(binned_image%n_nu,binned_image%n_ap,binned_image%n_view,binned_image%n_orig,binned_image%n_stokes))

          call mpi_reduce(binned_image%sed, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
          binned_image%sed = cube5d

          if(binned_image%uncertainties) then

             call mpi_reduce(binned_image%sed2, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             binned_image%sed2 = cube5d

             call mpi_reduce(binned_image%sedn, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             binned_image%sedn = cube5d

          end if

          deallocate(cube5d)

       end if

       if(binned_image%compute_image) then

          allocate(cube6d(binned_image%n_nu,binned_image%n_x,binned_image%n_y,binned_image%n_view,binned_image%n_orig,binned_image%n_stokes))

          call mpi_reduce(binned_image%img, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
          binned_image%img = cube6d

          if(binned_image%uncertainties) then

             call mpi_reduce(binned_image%img2, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             binned_image%img2 = cube6d

             call mpi_reduce(binned_image%imgn, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             binned_image%imgn = cube6d

          end if

          deallocate(cube6d)

       end if

    end if

    if(make_peeled_images) then

       do ip=1,n_groups

          if(peeled_image(ip)%compute_sed) then

             allocate(cube5d(peeled_image(ip)%n_nu,peeled_image(ip)%n_ap,peeled_image(ip)%n_view,peeled_image(ip)%n_orig, peeled_image(ip)%n_stokes))

             call mpi_reduce(peeled_image(ip)%sed, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             peeled_image(ip)%sed = cube5d

             if(peeled_image(ip)%uncertainties) then

                call mpi_reduce(peeled_image(ip)%sed2, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
                peeled_image(ip)%sed2 = cube5d

                call mpi_reduce(peeled_image(ip)%sedn, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
                peeled_image(ip)%sedn = cube5d

             end if

             deallocate(cube5d)

          end if

          if(peeled_image(ip)%compute_image) then

             allocate(cube6d(peeled_image(ip)%n_nu,peeled_image(ip)%n_x,peeled_image(ip)%n_y,peeled_image(ip)%n_view,peeled_image(ip)%n_orig,peeled_image(ip)%n_stokes))

             call mpi_reduce(peeled_image(ip)%img, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
             peeled_image(ip)%img = cube6d

             if(peeled_image(ip)%uncertainties) then

                call mpi_reduce(peeled_image(ip)%img2, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
                peeled_image(ip)%img2 = cube6d

                call mpi_reduce(peeled_image(ip)%imgn, cube6d, size(cube6d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
                peeled_image(ip)%imgn = cube6d

             end if

             deallocate(cube6d)

          end if

       end do

    end if

  end subroutine mp_collect_images

end module mpi_routines
