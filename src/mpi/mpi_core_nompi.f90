module mpi_core

  implicit none
  save

  private
  public :: mp_initialize
  public :: mp_stop
  public :: mp_join
  public :: main_process

  ! Whether MPI is enabled
  logical :: mpi_enabled = .false.

  ! Rank of current process
  integer, parameter, public :: rank = 0

contains

  logical function main_process()
    main_process = .true.
  end function main_process

  subroutine mp_initialize()
  end subroutine mp_initialize

  subroutine mp_stop()
  end subroutine mp_stop

  subroutine mp_join()
  end subroutine mp_join

end module mpi_core
