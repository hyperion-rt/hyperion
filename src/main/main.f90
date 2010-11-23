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

  call hdf5_set_compression(.true.)

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

  if(main_process()) call message_section('Initial setup')

  handle_in = hdf5_open_read(input_file)

  ! Prepare output directory
  if(main_process()) then
     handle_out = hdf5_open_new(output_file)
  end if

  ! Wait for all threads
  call mp_join()

  call cpu_time(time1)

  call set_verbose_level(100)
  call mp_set_random_seed()
  call setup_initial(handle_in)

  ! Wait for all threads
  call mp_join()  

  if(n_dust==0 .and. n_lucy_iter > 0) then
     call warn("main", "no dust present, so skipping temperature iterations")
     n_lucy_iter=0
  end if

  ! Loop over Lucy iterations
  do iter=1,n_lucy_iter

     ! Display message
     if(main_process()) call message_number(1,' [main] starting Lucy iteration ',iter,'(I0)','')

     ! Do the RT
     call do_lucy(n_lucy_photons, n_stats)

     ! Wait for all threads
     call mp_join()

     ! Display message
     if(main_process()) call message(1,' [main] exiting Lucy iteration')

     ! Output files
     if(main_process()) call output_grid(handle_out, iter, n_lucy_iter)

  end do

  ! Set up image-related variables
  call setup_final_iteration(handle_in)

  ! Prepare output directories for images/SEDs
  if(main_process()) then
     if(make_binned_images) g_binned = hdf5_create_group(handle_out, 'Binned')
     if(make_peeled_images) g_peeled = hdf5_create_group(handle_out, 'Peeled')
  end if

  ! Display message
  if(main_process()) call message(1,' [main] starting final iteration')

  ! Do the RT
  if(use_exact_nu) then
     call do_final_mono(n_last_photons, 20*n_last_photons, n_stats, use_raytracing)
  else
     call do_final(n_last_photons, n_stats, use_raytracing)
  end if

  ! Display message
  if(main_process()) call message(1,' [main] exiting final iteration')

  if(use_raytracing) then

     ! Display message
     if(main_process()) call message(1,' [main] starting raytracing iteration')

     ! Do the raytracing
     call do_raytracing(n_raytracing_photons_star,n_raytracing_photons_dust, n_stats)

     ! Display message
     if(main_process()) call message(1,' [main] exiting raytracing iteration')

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
     write(*,'(" Total CPU time elapsed: ",F16.2)') time
  end if

  ! Stop multi-processing
  call mp_stop()

end program main
