! Explicit interfaces for the MPI routines that take "choice" buffer
! arguments (buffers that can be of any type and rank). The mpich "mpi"
! module does not declare interfaces for these routines, so gfortran 10 and
! later warns about every pair of calls to the same routine with buffers of
! different types or ranks. The NO_ARG_CHECK attribute (supported by both
! gfortran and the Intel compilers, and ignored as a comment by others)
! reproduces the "ignore type, kind and rank" behaviour that MPI
! implementations themselves use for choice buffers, while keeping the
! non-buffer arguments type-checked.

module mpi_interfaces

  implicit none

  interface

     subroutine mpi_bcast(buffer, count, datatype, root, comm, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: buffer
       type(*) :: buffer(*)
       integer :: count, datatype, root, comm, ierror
     end subroutine mpi_bcast

     subroutine mpi_reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
       type(*) :: sendbuf(*), recvbuf(*)
       integer :: count, datatype, op, root, comm, ierror
     end subroutine mpi_reduce

     subroutine mpi_allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
       type(*) :: sendbuf(*), recvbuf(*)
       integer :: count, datatype, op, comm, ierror
     end subroutine mpi_allreduce

     subroutine mpi_recv(buf, count, datatype, source, tag, comm, status, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: buf
       type(*) :: buf(*)
       integer :: count, datatype, source, tag, comm, status(*), ierror
     end subroutine mpi_recv

     subroutine mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: buf
       type(*) :: buf(*)
       integer :: count, datatype, dest, tag, comm, request, ierror
     end subroutine mpi_isend

     subroutine mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror)
       !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: buf
       type(*) :: buf(*)
       integer :: count, datatype, source, tag, comm, request, ierror
     end subroutine mpi_irecv

  end interface

end module mpi_interfaces
