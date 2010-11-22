module performance

  implicit none
  save

  private
  public :: perf_init
  public :: perf_display

  real :: time_1,time_2
  ! time interval to compute the statistics

contains

  subroutine perf_init
    implicit none
    call cpu_time(time_1)
    write(*,*)
    write(*,'("   # Photons    CPU time (sec)    Photons/sec  ")')
    write(*,'(" ----------------------------------------------")')
  end subroutine perf_init

  subroutine perf_display(count)
    implicit none
    integer,intent(in) :: count
    call cpu_time(time_2)
    write(*,'(1X,3X,I7,3X,4X,F10.1,4X,4X,F9.2,4X)') count,time_2-time_1,real(count)/(time_2-time_1+1.e-30)
  end subroutine perf_display

end module performance