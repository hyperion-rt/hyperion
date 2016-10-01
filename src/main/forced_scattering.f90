module forced_scattering

  use core_lib, only : dp, random

  implicit none
  save

  private
  public :: forced_scattering_wr99
  public :: forced_scattering_baes16

  ! The optical depth at which we transition to approximating (1 - exp(-tau))
  ! by simply tau
  real(dp), parameter :: TAU_THRES = 1.e-7_dp

contains

  subroutine forced_scattering_wr99(tau_escape, tau, weight)

    ! Simple forced first scattering from Wood & Reynolds, 1999, The
    ! Astrophysical Journal, 525, 799:
    !
    ! http://dx.doi.org/10.1086/307939
    !
    ! This changes the probability density function for tau from exp(-tau)
    ! defined between 0 and infinity to being a truncated decaying exponential
    ! going from 0 to tau_escape. The PDF is thus:
    !
    !
    ! PDF = exp(-tau) / (1 - exp(-tau_escape))
    !
    ! and the CDF is
    !
    ! CDF = (1 - exp(-tau)) / (1 - exp(-tau_escape))

    implicit none

    real(dp),intent(in) :: tau_escape
    real(dp),intent(out) :: tau, weight
    real(dp) :: xi, one_minus_exp

    call random(xi)

    if(tau_escape > TAU_THRES) then
      one_minus_exp = (1._dp - exp(-tau_escape))
    else
      one_minus_exp = tau_escape
    end if

    tau = -log(1._dp - xi * one_minus_exp)
    weight = one_minus_exp

  end subroutine forced_scattering_wr99

  subroutine forced_scattering_baes16(eta, tau_escape, tau, weight)

    ! Forced first scattering with composite biasing from Baes et al. 2019,
    ! Astronomy and Astrophysics, 590, A55:
    !
    ! http://dx.doi.org/10.1086/307939
    !
    ! This changes the probability density function for tau from exp(-tau)
    ! defined between 0 and infinity to being a composite of a truncated
    ! decaying exponential and a constant. The PDF is thus:
    !
    !
    ! PDF = (1 - eta) * exp(-tau) / (1 - exp(-tau_escape)) + eta / tau_escape
    !
    ! where eta is in the range [0:1], and the CDF is
    !
    ! CDF = (1 - eta) * (1 - exp(-tau)) / (1 - exp(-tau_escape)) + eta * tau / tau_escape
    !
    ! For simplicity, we define:
    !
    ! alpha = (1 - eta) / (1 - exp(-tau_escape))
    ! beta = eta / tau_escape
    !
    ! so that the CDF becomes:
    !
    ! CDF = alpha * (1 - exp(-tau)) + beta * tau
    !
    ! Unfortunately, we cannot solve CDF = xi analytically for tau (where xi is
    ! a random number in the range [0:1]), so instead we solve this by searching
    ! the CDF by bisection.

    implicit none

    real(dp),intent(in) :: eta, tau_escape
    real(dp),intent(out) :: tau, weight
    real(dp) :: alpha, beta, tau_min, tau_max, xi, xi_test, one_minus_exp
    integer :: i

    if(tau_escape > TAU_THRES) then
      one_minus_exp = 1._dp - exp(-tau_escape)
    else
      one_minus_exp = tau_escape
    end if

    ! Pre-define alpha and beta for convenience and performance
    alpha = (1._dp - eta) / one_minus_exp
    beta = eta / tau_escape

    ! Search tau by bisection - we know in advance it will be in the range
    ! 0 to tau_escape, so we use this as initial limits
    tau_min = 0._dp
    tau_max = tau_escape

    call random(xi)

    ! 50 iterations is enough to get a result to machine precision
    do i=1,60
       tau = 0.5_dp * (tau_min + tau_max)
       if (tau > TAU_THRES) then
         xi_test = alpha * (1._dp - exp(-tau)) + beta * tau
       else
         xi_test = alpha * tau + beta * tau
       end if
       if (xi_test > xi) then
          tau_max = tau
       else
          tau_min = tau
       end if
    end do

    tau = 0.5_dp * (tau_min + tau_max)

    weight = 1._dp / (alpha + beta * exp(tau))

  end subroutine forced_scattering_baes16

end module forced_scattering
