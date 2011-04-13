program main

  use mpi_core
  use mpi_routines
  use setup
  use binned_images
  use peeled_images
  use iteration_lucy
  use iteration_final
  use iteration_final_mono
  use iteration_raytracing
  use grid_generic
  use settings

  implicit none

  integer :: iter
  character(len=1000) :: input_file, output_file
  integer(hid_t) :: handle_in, handle_out, g_peeled, g_binned
  real(dp) :: time1, time2, time
  logical :: converged
  character(len=30) :: datetime

  character(len=5), parameter :: version = '0.7.6'

  call mp_set_compression(.true.)

  ! Retrieve command-line arguments
  call get_command_argument(1, input_file)
  call get_command_argument(2, output_file)

  ! Check that both arguments were given
  if(trim(input_file)=="".or.trim(output_file)=="") then
     stop "Usage: bin/rt input_file output_file"
  end if

  ! Start up multi-processing if needed
  call mp_initialize()

  ! SETUP

  call check_file_exists(input_file)
  handle_in = mp_open_read(input_file)

  ! Prepare output directory
  if(main_process()) then
     handle_out = mp_open_new(output_file)
     call mp_write_keyword(handle_out, '/', 'fortran_version', version)
  end if

  if(main_process()) then
     write(*,*) repeat('-',60)
     datetime = now()
     write(*,'(" Started on ",A)') trim(datetime)
     call mp_write_keyword(handle_out, '/', 'date_started', trim(datetime))
     write(*,'(" Input:  ", A)') trim(input_file)
     write(*,'(" Output: ", A)') trim(output_file)
     write(*,*) repeat('-',60)
  end if

  ! Wait for all threads
  call mp_join()

  call cpu_time(time1)

  call set_verbose_level(100)
  call mp_set_random_seed()
  call setup_initial(handle_in)

  ! Wait for all threads
  call mp_join()

  ! Loop over Lucy iterations
  do iter=1,n_lucy_iter

     ! Display message
     if(main_process()) write(*,'(" [main] starting Lucy iteration ", I0)') iter

     ! Do the RT
     call do_lucy(n_lucy_photons, n_stats)

     ! Wait for all threads
     call mp_join()

     ! Display message
     if(main_process()) write(*,'(" [main] exiting Lucy iteration")')

     ! Check for convergence
     if(check_convergence) then

        if(main_process()) converged = specific_energy_abs_converged()

        call mp_broadcast_convergence(converged)

        if(converged) then

           if(main_process()) then

              write(*,'("      ------ Temperature calculation converged -----")')
              write(*,*)

              ! Output files (and signal that this is the last iteration)
              call output_grid(handle_out, iter, iter)

           end if

           ! Exit the temperature iteration
           exit

        end if

     end if

     ! Output files. The following needs to be executed on all ranks because
     ! the MPI AMR version needs to sync during mp_path_exists.
     call output_grid(handle_out, iter, n_lucy_iter)

  end do

  if(main_process()) then
     call mp_write_keyword(handle_out, '/', 'converged', converged)
     if(converged) then
        call mp_write_keyword(handle_out, '/', 'iterations', iter)
     else
        call mp_write_keyword(handle_out, '/', 'iterations', n_lucy_iter)
     end if
  end if

  ! Set up image-related variables
  call setup_final_iteration(handle_in)

  ! Prepare output directories for images/SEDs
  if(main_process()) then
     if(make_binned_images) g_binned = mp_create_group(handle_out, 'Binned')
     if(make_peeled_images) g_peeled = mp_create_group(handle_out, 'Peeled')
  end if

  ! Display message
  if(main_process()) write(*,'(" [main] starting final iteration")')

  ! Do the RT
  if(use_exact_nu) then
     call do_final_mono(n_last_photons_sources, n_last_photons_dust, n_stats, use_raytracing)
  else
     call do_final(n_last_photons, n_stats, use_raytracing)
  end if

  ! Display message
  if(main_process()) write(*,'(" [main] exiting final iteration")')

  if(use_raytracing) then

     ! Display message
     if(main_process()) write(*,'(" [main] starting raytracing iteration")')

     ! Do the raytracing
     call do_raytracing(n_raytracing_photons_sources,n_raytracing_photons_dust, n_stats)

     ! Display message
     if(main_process()) write(*,'(" [main] exiting raytracing iteration")')

  end if

  call mp_collect_results()

  ! Write out images
  if(main_process()) then
     if(make_binned_images) call binned_images_write(g_binned)
     if(make_peeled_images) call peeled_images_write(g_peeled)
  end if

  call cpu_time(time2)

  time = time2 - time1

  call mp_sync_cputime(time)

  call mp_join()

  if(main_process()) then
     write(*,*) repeat('-',60)
     write(*,'(" Total CPU time elapsed: ",F16.2)') time
     call mp_write_keyword(handle_out, '/', 'cpu_time', time)
     datetime = now()
     write(*,'(" Ended on ",A)') trim(datetime)
     call mp_write_keyword(handle_out, '/', 'date_ended', trim(datetime))
     write(*,*) repeat('-',60)
  end if

  ! Close output file
  call mp_close(handle_out)

  ! Stop multi-processing
  call mp_stop()

end program main
