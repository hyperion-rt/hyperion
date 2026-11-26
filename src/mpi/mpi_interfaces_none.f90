! Empty version of the mpi_interfaces module, used instead of
! mpi_interfaces.f90 when the "mpi" module of the MPI library already declares
! interfaces for the MPI routines that take choice buffers (as Open MPI does),
! since declaring them again makes the names ambiguous. The configure script
! chooses between the two.

module mpi_interfaces

  implicit none

end module mpi_interfaces
