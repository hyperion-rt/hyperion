program test

  use mpi

  implicit none

  integer :: ierr, rank, nproc

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

  print *, 'MPI Installation Successful!'

end program test

