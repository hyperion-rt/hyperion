module mpi_core

  use iso_c_binding, only : c_char, c_int, c_size_t, c_null_char

  implicit none
  save

  private
  public :: mp_initialize
  public :: mp_stop
  public :: mp_join
  public :: main_process
  public :: mp_node_name

  ! Whether MPI is enabled
  logical, public :: mpi_enabled = .false.

  ! Rank of current process
  integer, parameter, public :: rank = 0

  interface
     integer(c_int) function gethostname(name, length) bind(c, name='gethostname')
       import :: c_char, c_int, c_size_t
       character(kind=c_char), dimension(*) :: name
       integer(c_size_t), value :: length
     end function gethostname
  end interface

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

  function mp_node_name() result(name)
    ! Return the name of the machine the code is running on
    implicit none
    character(len=256) :: name
    character(kind=c_char), dimension(257) :: buffer
    integer :: i
    name = ''
    buffer = c_null_char
    if(gethostname(buffer, int(size(buffer) - 1, c_size_t)) == 0) then
       do i = 1, len(name)
          if(buffer(i) == c_null_char) exit
          name(i:i) = buffer(i)
       end do
    end if
    if(len_trim(name) == 0) name = 'unknown'
  end function mp_node_name

end module mpi_core
