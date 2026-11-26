module mpi_core

  use mpi

  implicit none
  save

  private
  public :: mp_initialize
  public :: mp_stop
  public :: mp_join
  public :: main_process
  public :: mp_node_name

  ! MPI errors
  integer :: ierr

  ! Whether MPI is enabled
  logical,public :: mpi_enabled = .true.

  ! Number of active processes
  integer,public :: nproc

  ! Rank of current process
  integer,public :: rank

  ! The rank of the main process
  integer,parameter,public :: rank_main = 0

contains

  logical function main_process()
    main_process = rank == rank_main
  end function main_process

  subroutine mp_initialize()
    implicit none
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
  end subroutine mp_initialize

  subroutine mp_stop()
    implicit none
    call mpi_finalize(ierr)
  end subroutine mp_stop

  subroutine mp_join()
    implicit none
    call mpi_barrier(mpi_comm_world, ierr)
  end subroutine mp_join

  function mp_node_name() result(name)
    ! Return the name of the node the current process is running on
    implicit none
    character(len=mpi_max_processor_name) :: name
    integer :: length
    name = ''
    call mpi_get_processor_name(name, length, ierr)
    if(len_trim(name) == 0) name = 'unknown'
  end function mp_node_name

end module mpi_core
