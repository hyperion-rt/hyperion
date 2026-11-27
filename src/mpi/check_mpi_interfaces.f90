! Program compiled by the configure script (but not used otherwise) to check
! whether the explicit interfaces in mpi_interfaces.f90 can be used together
! with the "mpi" module of the MPI library. If the "mpi" module already
! declares these interfaces (as Open MPI does), the names are ambiguous and
! this does not compile, in which case mpi_interfaces_none.f90 is used instead.

program check_mpi_interfaces

  use mpi
  use mpi_interfaces

  implicit none

  integer :: a(1), ierr

  call mpi_bcast(a, 1, mpi_integer, 0, mpi_comm_world, ierr)

end program check_mpi_interfaces
