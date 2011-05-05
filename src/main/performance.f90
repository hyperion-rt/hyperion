module performance

  use core_lib

  implicit none
  save

  private
  public :: perf_header
  public :: perf_footer
  public :: perf_numbers

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
    if(time < tiny(1._dp)) then
       write(*,'(1X,3X,I12,3X,4X,F10.1,4X,4X,"   ...   ")') count,time
    else
       write(*,'(1X,3X,I12,3X,4X,F10.1,4X,4X,F9.2,4X)') count,time,real(count)/(time)
    end if
  end subroutine perf_numbers

end module performance
