module mpi_routines

  use mpi
  use mpi_core
  use core_lib
  use grid_physics
  use sources, only : energy_current
  use setup

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
  public :: mp_collect
  public :: mp_broadcast_temperature
  public :: mp_collect_results
  public :: mp_sync_energy
  public :: mp_sync_cputime
  public :: mp_set_random_seed

  real(dp) :: time_curr
  integer(idp) :: n_completed

contains

  subroutine mp_reset_first()
    implicit none
    first = .true.
    time_curr = 0._dp
    n_completed = 0
  end subroutine mp_reset_first

  subroutine mp_n_photons(n_photons_tot, n_photons_curr, n_photons_chunk, n_photons)

    implicit none

    ! Total number of photons, and size of a chunk
    integer(idp),intent(in) :: n_photons_tot, n_photons_chunk

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
             if(n_photons_tot - n_photons_curr <= n_photons_chunk*2_idp) then
                n_photons_send(ir) = n_photons_chunk / 10_idp
             else
                n_photons_send(ir) = n_photons_chunk
             end if
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
          n_photons = n_photons_chunk / 10_idp
          n_photons_curr = n_photons_curr + n_photons

       end if

       if(first) first=.false.       

       if(mod(n_photons_curr,n_photons_chunk)==0) then
          write(*,'(1X,3X,I12,3X,4X,F10.1,4X,4X,F9.2,4X)') n_photons_curr,time_curr,dble(n_photons_curr)/time_curr
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

  subroutine mp_set_random_seed()
    implicit none
    call set_seed(-124902+rank)
  end subroutine mp_set_random_seed

  subroutine mp_collect()

    implicit none
    real(dp) :: tmp
    real(dp) :: dummy_dp
    integer(idp) :: dummy_idp
    real(dp),allocatable :: tmp_2d(:,:)
    integer(idp),allocatable :: tmp_int_1d(:)

    call mpi_reduce(energy_current, tmp, 1, mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
    energy_current = tmp

    if(main_process()) then
       allocate(tmp_2d(size(specific_energy_abs,1),size(specific_energy_abs,2)))
       call mpi_reduce(specific_energy_abs, tmp_2d, size(specific_energy_abs), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
       specific_energy_abs = tmp_2d
       deallocate(tmp_2d)
    else
       call mpi_reduce(specific_energy_abs, dummy_dp, size(specific_energy_abs), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr)
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

  end subroutine mp_collect

  subroutine mp_broadcast_temperature() 

    implicit none

    call mpi_bcast(temperature, size(temperature), mpi_real8, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(specific_energy_abs, size(specific_energy_abs), mpi_real8, rank_main, mpi_comm_world, ierr)
    call mpi_bcast(energy_abs_tot, size(energy_abs_tot), mpi_real8, rank_main, mpi_comm_world, ierr)

    if(allocated(temperature_mean)) then
       call mpi_bcast(temperature_mean, size(temperature_mean), mpi_real8, rank_main, mpi_comm_world, ierr)
    end if

  end subroutine mp_broadcast_temperature

  subroutine mp_sync_energy()
    implicit none
    real(dp) :: tmp
    call mpi_allreduce(energy_current, tmp, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    energy_current = tmp
  end subroutine mp_sync_energy

  subroutine mp_sync_cputime(cputime)
    implicit none
    real(dp),intent(inout) :: cputime
    real(dp) :: tmp
    call mpi_allreduce(cputime, tmp, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    cputime = tmp
  end subroutine mp_sync_cputime

  subroutine mp_collect_results()

    use binned_images
    use peeled_images

    implicit none

    integer :: ip

    real(dp),allocatable :: cube4d(:,:,:,:)
    real(dp),allocatable :: cube5d(:,:,:,:,:)

    call mp_sync_energy()

    if(make_binned_images) then

       allocate(cube4d(binned_image%n_nu,binned_image%n_ap,binned_image%n_view,binned_image%n_orig))
       allocate(cube5d(binned_image%n_x,binned_image%n_y,binned_image%n_nu,binned_image%n_view,binned_image%n_orig))

       call mpi_reduce(binned_image%sed%i, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed%i = cube4d
       call mpi_reduce(binned_image%sed%q, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed%q = cube4d
       call mpi_reduce(binned_image%sed%u, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed%u = cube4d
       call mpi_reduce(binned_image%sed%v, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed%v = cube4d

       if(binned_image%uncertainties) then

          call mpi_reduce(binned_image%sed2%i, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed2%i = cube4d
          call mpi_reduce(binned_image%sed2%q, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed2%q = cube4d
          call mpi_reduce(binned_image%sed2%u, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed2%u = cube4d
          call mpi_reduce(binned_image%sed2%v, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sed2%v = cube4d

          call mpi_reduce(binned_image%sedn, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%sedn = cube4d

       end if

       call mpi_reduce(binned_image%img%i, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img%i = cube5d
       call mpi_reduce(binned_image%img%q, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img%q = cube5d
       call mpi_reduce(binned_image%img%u, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img%u = cube5d
       call mpi_reduce(binned_image%img%v, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img%v = cube5d

       if(binned_image%uncertainties) then

          call mpi_reduce(binned_image%img2%i, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img2%i = cube5d
          call mpi_reduce(binned_image%img2%q, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img2%q = cube5d
          call mpi_reduce(binned_image%img2%u, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img2%u = cube5d
          call mpi_reduce(binned_image%img2%v, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%img2%v = cube5d

          call mpi_reduce(binned_image%imgn, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; binned_image%imgn = cube5d

       end if

       deallocate(cube4d)
       deallocate(cube5d)

    end if

    if(make_peeled_images) then

       do ip=1,n_groups

          allocate(cube4d(peeled_image(ip)%n_nu,peeled_image(ip)%n_ap,peeled_image(ip)%n_view,peeled_image(ip)%n_orig))
          allocate(cube5d(peeled_image(ip)%n_x,peeled_image(ip)%n_y,peeled_image(ip)%n_nu,peeled_image(ip)%n_view,peeled_image(ip)%n_orig))

          call mpi_reduce(peeled_image(ip)%sed%i, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed%i = cube4d
          call mpi_reduce(peeled_image(ip)%sed%q, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed%q = cube4d
          call mpi_reduce(peeled_image(ip)%sed%u, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed%u = cube4d
          call mpi_reduce(peeled_image(ip)%sed%v, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed%v = cube4d

          if(peeled_image(ip)%uncertainties) then

             call mpi_reduce(peeled_image(ip)%sed2%i, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed2%i = cube4d
             call mpi_reduce(peeled_image(ip)%sed2%q, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed2%q = cube4d
             call mpi_reduce(peeled_image(ip)%sed2%u, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed2%u = cube4d
             call mpi_reduce(peeled_image(ip)%sed2%v, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sed2%v = cube4d

             call mpi_reduce(peeled_image(ip)%sedn, cube4d, size(cube4d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%sedn = cube4d

          end if

          call mpi_reduce(peeled_image(ip)%img%i, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img%i = cube5d
          call mpi_reduce(peeled_image(ip)%img%q, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img%q = cube5d
          call mpi_reduce(peeled_image(ip)%img%u, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img%u = cube5d
          call mpi_reduce(peeled_image(ip)%img%v, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img%v = cube5d

          if(peeled_image(ip)%uncertainties) then

             call mpi_reduce(peeled_image(ip)%img2%i, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img2%i = cube5d
             call mpi_reduce(peeled_image(ip)%img2%q, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img2%q = cube5d
             call mpi_reduce(peeled_image(ip)%img2%u, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img2%u = cube5d
             call mpi_reduce(peeled_image(ip)%img2%v, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%img2%v = cube5d

             call mpi_reduce(peeled_image(ip)%imgn, cube5d, size(cube5d), mpi_real8, mpi_sum, rank_main, mpi_comm_world, ierr) ; peeled_image(ip)%imgn = cube5d

          end if

          deallocate(cube4d)
          deallocate(cube5d)

       end do

    end if

  end subroutine mp_collect_results

end module mpi_routines
