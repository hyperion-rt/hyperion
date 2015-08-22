program main

  use core_lib, only : hid_t, dp, now, check_file_exists, warn, set_verbose_level

  use mpi_core, only : main_process, mp_join, mp_stop, mp_initialize

  use mpi_routines, only : mp_sync, &
       &                   mp_set_random_seed, &
       &                   mp_broadcast_convergence, &
       &                   mp_collect_images

  use mpi_hdf5_io, only : mp_open_read, &
       &                  mp_open_new, &
       &                  mp_test_version, &
       &                  mp_create_group, &
       &                  mp_set_compression, &
       &                  mp_write_keyword, &
       &                  mp_exists_keyword, &
       &                  mp_read_keyword, &
       &                  mp_copy_group, &
       &                  mp_create_external_link, &
       &                  mp_close, &
       &                  mp_finalize

  use setup, only : setup_initial, &
       &            setup_final_iteration

  use binned_images, only : make_binned_images, &
       &                    binned_images_write
  use peeled_images, only : make_peeled_images, &
       &                    peeled_images_write

  use iteration_lucy, only : do_lucy
  use iteration_final, only : do_final
  use iteration_final_mono, only : do_final_mono
  use iteration_raytracing, only : do_raytracing

  use grid_generic, only : output_grid
  use grid_physics, only : specific_energy_converged

  use settings, only : n_initial_iter, &
       &               n_initial_photons, &
       &               n_last_photons_sources, &
       &               n_last_photons_dust, &
       &               n_last_photons, &
       &               n_raytracing_photons_sources, &
       &               n_raytracing_photons_dust, &
       &               n_stats, &
       &               check_convergence, &
       &               use_exact_nu, &
       &               use_raytracing

  use counters, only : killed_photons_geo, &
       &               killed_photons_int

  implicit none

  integer :: iter
  character(len=1000) :: command, force_option, input_file, output_file
  integer(hid_t) :: handle_in, handle_out, g_peeled, g_binned, g_input
  real(dp) :: time1, time2, time
  logical :: converged, copy_input
  character(len=30) :: datetime
  integer :: seed
  character(len=100) :: group_name
  integer(hid_t) :: handle_iter
  logical :: confirm_overwrite, show_usage
  integer :: n_args

  character(len=5), parameter :: fortran_version = '0.9.8'

  call mp_set_compression(.true.)

  n_args = command_argument_count()

  ! Retrieve command-line arguments
  confirm_overwrite = .true.
  show_usage = .false.

  call get_command_argument(0, command)
  if(n_args == 2) then
     call get_command_argument(1, input_file)
     call get_command_argument(2, output_file)
  else if(n_args == 3) then
     call get_command_argument(1, force_option)
     call get_command_argument(2, input_file)
     call get_command_argument(3, output_file)
     if(trim(force_option) == '-f') then
        confirm_overwrite = .false.
     else
        show_usage = .true.
     end if
  else
     show_usage = .true.
  end if


  ! Start up multi-processing if needed
  call mp_initialize()

  ! Check that both arguments were given
  if(show_usage) then
     if(main_process()) write(*,'("Usage: ", A, " [-f] input_file output_file")') trim(command)
     call mp_stop()
     stop
  end if

  ! SETUP

  if(main_process()) then
     write(*,*) repeat('-',60)
     datetime = now()
     write(*,'(" Hyperion v",A)') fortran_version
     write(*,'(" Started on ",A)') trim(datetime)
     write(*,'(" Input:  ", A)') trim(input_file)
     write(*,'(" Output: ", A)') trim(output_file)
     write(*,*) repeat('-',60)
  end if

  if(main_process()) call check_file_exists(input_file)

  ! Wait for all threads
  call mp_join()

  handle_in = mp_open_read(input_file)

  ! Prepare output directory
  if(main_process()) then
     handle_out = mp_open_new(output_file, confirm=confirm_overwrite)
     call mp_write_keyword(handle_out, '/', 'date_started', trim(datetime))
     call mp_write_keyword(handle_out, '/', 'fortran_version', fortran_version)
  end if

  ! Include the input in the output file
  call mp_read_keyword(handle_in, '/', 'copy_input', copy_input)

  ! Add HDF5 version check - copy requires a recent version of HDF5 because
  ! writing tables with uneven field names (which we need) was buggy before
  ! that version.
  if(copy_input.and..not.mp_test_version(1, 8, 6)) then
     if(main_process()) call warn("main","copy_input option requires HDF5 1.8.6 or later, linking input")
     copy_input = .false.
  end if

  if(copy_input) then
     g_input = mp_create_group(handle_out, '/Input')
     call mp_copy_group(handle_in, '/', g_input, '.')
  else
     call mp_create_external_link(handle_out, '/Input', input_file, '/')
  end if

  ! Wait for all threads
  call mp_join()

  call cpu_time(time1)

  if(mp_exists_keyword(handle_in, '/', 'seed')) then
     call mp_read_keyword(handle_in, '/', 'seed', seed)
  else
     seed = -124902  ! value used before customized seed was implemented
  end if

  if(main_process()) write(*, '(" [main] using random seed = ", I0)') seed

  call set_verbose_level(100)
  call mp_set_random_seed(seed)
  call setup_initial(handle_in)

  ! Wait for all threads
  call mp_join()

  ! Loop over Lucy iterations
  do iter=1,n_initial_iter

     ! Reset killed photon counters to zero
     killed_photons_geo = 0
     killed_photons_int = 0

     ! Display message
     if(main_process()) write(*,'(" [main] starting Lucy iteration ", I0)') iter

     ! Do the RT
     call do_lucy(n_initial_photons, n_stats)

     ! Wait for all threads
     call mp_join()

     ! Display message
     if(main_process()) write(*,'(" [main] exiting Lucy iteration")')

     ! Check for convergence
     if(check_convergence) then

        if(main_process()) converged = specific_energy_converged()

        call mp_broadcast_convergence(converged)

        if(converged) then

           if(main_process()) then
              write(*,'("      ------ Specific energy calculation converged -----")')
              write(*,*)
           end if

        end if

     end if

     ! Create HDF5 group name
     write(group_name, '("iteration_",I5.5)') iter
     handle_iter = mp_create_group(handle_out, group_name)

     ! Output files. The following needs to be executed on all ranks because
     ! the MPI AMR version needs to sync during mp_path_exists.
     if(check_convergence.and.converged) then
        call output_grid(handle_iter, iter, iter)
     else
        call output_grid(handle_iter, iter, n_initial_iter)
     end if

     ! Sync killed photon counters
     call mp_sync(killed_photons_geo)
     call mp_sync(killed_photons_int)

     ! Write out killed photon information
     call mp_write_keyword(handle_iter, '.', 'killed_photons_geo', killed_photons_geo)
     call mp_write_keyword(handle_iter, '.', 'killed_photons_int', killed_photons_int)

     if(check_convergence.and.converged) exit

  end do

  ! CONVERGENCE INFORMATION

  if(main_process()) then

     ! Write out convergence information
     call mp_write_keyword(handle_out, '/', 'converged', converged)
     if(converged) then
        call mp_write_keyword(handle_out, '/', 'iterations', iter)
     else
        call mp_write_keyword(handle_out, '/', 'iterations', n_initial_iter)
     end if

  end if

  ! KILLED PHOTON INFORMATION

  ! Reset killed photon counters to zero
  killed_photons_geo = 0
  killed_photons_int = 0

  ! FINAL ITERATION

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

  ! KILLED PHOTON INFORMATION

  ! Sync killed photon counters
  call mp_sync(killed_photons_geo)
  call mp_sync(killed_photons_int)

  ! Write out killed photon information
  call mp_write_keyword(handle_out, '/', 'killed_photons_geo_final', killed_photons_geo)
  call mp_write_keyword(handle_out, '/', 'killed_photons_int_final', killed_photons_int)

  ! Reset killed photon counters to zero
  killed_photons_geo = 0
  killed_photons_int = 0

  ! RAYTRACING ITERATION

  if(use_raytracing) then

     ! Display message
     if(main_process()) write(*,'(" [main] starting raytracing iteration")')

     ! Do the raytracing
     call do_raytracing(n_raytracing_photons_sources,n_raytracing_photons_dust, n_stats)

     ! Display message
     if(main_process()) write(*,'(" [main] exiting raytracing iteration")')

  end if

  ! KILLED PHOTON INFORMATION

  ! Sync killed photon counters
  call mp_sync(killed_photons_geo)
  call mp_sync(killed_photons_int)

  ! Write out killed photon information
  call mp_write_keyword(handle_out, '/', 'killed_photons_geo_raytracing', killed_photons_geo)
  call mp_write_keyword(handle_out, '/', 'killed_photons_int_raytracing', killed_photons_int)

  ! OUTPUT

  ! Collect images (and SEDs) to the main process
  call mp_collect_images()

  ! Write out images
  if(main_process()) then
     if(make_binned_images) call binned_images_write(g_binned)
     if(make_peeled_images) call peeled_images_write(g_peeled)
  end if

  call cpu_time(time2)

  time = time2 - time1

  call mp_sync(time)

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

  ! Finalize HDF5
  call mp_finalize()

  ! Stop multi-processing
  call mp_stop()

end program main
