module performance

  use core_lib

  implicit none
  save

  private
  public :: perf_header
  public :: perf_footer
  public :: perf_numbers
  public :: perf_init
  public :: perf_display

  real(dp) :: time_1,time_2
  ! time interval to compute the statistics

contains

  subroutine perf_header()
    implicit none
    write(*,*)
    write(*,'("        # Photons    CPU time (sec)    Photons/sec  ")')
    write(*,'("      ----------------------------------------------")')
  end subroutine perf_header

  subroutine perf_footer()
    implicit none
    write(*,'("      ----------------------------------------------")')
    write(*,*)
  end subroutine perf_footer

  subroutine perf_numbers(count, time)
    implicit none
    integer(idp),intent(in) :: count
    real(dp),intent(in) :: time
    write(*,'(1X,3X,I12,3X,4X,F10.1,4X,4X,F9.2,4X)') count,time,real(count)/(time+1.e-30)
  end subroutine perf_numbers

  subroutine perf_init()
    implicit none
    call cpu_time(time_1)
    call perf_header()
  end subroutine perf_init

  subroutine perf_display(count)
    implicit none
    integer(idp),intent(in) :: count
    call cpu_time(time_2)
    call perf_numbers(count, time_2 - time_1)
  end subroutine perf_display

end module performance
