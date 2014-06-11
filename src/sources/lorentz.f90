module lorentz

  use base_types
  use type_vector3d
  use type_angle3d

  implicit none
  save

  private
  public :: doppler_shift

contains

  function doppler_shift(nu0, a, v) result(nu)

    ! Relativistic doppler shift without beaming - that is, the frequency
    ! changes but the direction of the photon remains the same. This is used
    ! because if the full Lorentz transform was used, the direction of
    ! propagation of the photons would change, which would make proper
    ! raytracing and peeloff, etc. much more difficult.
    !
    ! Parameters
    ! ----------
    ! nu0 : real(dp)
    !     The frequency of the photon in the original frame of reference.
    ! a : angle3d_dp
    !     The direction vector of the photon in the original frame of reference
    ! v : vector3d_dp
    !     The velocity of the frame of reference relative to the original frame
    !     of reference
    !
    ! Returns
    ! -------
    ! nu : real(dp)
    !     The new frequency for the photon

    real(dp), intent(in) :: nu0
    type(angle3d_dp),intent(in) :: a
    type(vector3d_dp),intent(in) :: v
    real(dp) :: nu

    type(vector3d_dp) :: d

    real(dp), parameter :: c = 29979245800.

    type(vector3d_dp) :: beta

    call angle3d_to_vector3d(a, d)

    beta = v / c
    nu = (1._dp + (beta .dot. d)) * nu0 / sqrt(1._dp - (beta .dot. beta))

  end function doppler_shift

  subroutine lorentz_boost(nu, a, v)

    ! Parameters
    ! ----------
    ! nu : real(dp)
    !     The frequency of the photon in the original frame of reference.
    ! a : angle3d_dp
    !     The direction vector of the photon in the original frame of reference
    ! v : vector3d_dp
    !     The velocity of the frame of reference relative to the original frame
    !     of reference
    !
    ! Returns
    ! -------
    ! nu : real(dp)
    !     The new frequency for the photon
    ! a : angle3d_dp
    !     The new direction vector for the photons

    real(dp), intent(inout) :: nu
    type(angle3d_dp),intent(inout) :: a
    type(vector3d_dp),intent(in) :: v

    type(vector3d_dp) :: d

    real(dp), parameter :: c = 29979245800.

    type(vector3d_dp) :: beta

    real(dp) :: gamma, beta_sq
    real(dp) :: L(4, 4), F1(4), F2(4)
    real(dp) :: norm

    ! TODO - optimize this because we are wasting time converting a to v and vice-versa
    call angle3d_to_vector3d(a, d)

    ! Compute relativistic quantities
    beta = v / c
    beta_sq = beta .dot. beta
    gamma = 1._8 - beta_sq

    ! Compute Lorentz transformation matrix

    L(1, 1) = gamma
    L(1, 2) = - gamma * beta%x
    L(1, 3) = - gamma * beta%y
    L(1, 4) = - gamma * beta%z

    L(2, 1) = L(1, 2)
    L(2, 2) = 1._8 + (gamma - 1) * beta%x * beta%x / beta_sq
    L(2, 3) = (gamma - 1) * beta%x * beta%y / beta_sq
    L(2, 4) = (gamma - 1) * beta%x * beta%z / beta_sq

    L(3, 1) = L(1, 3)
    L(3, 2) = L(2, 3)
    L(3, 3) = 1._8 + (gamma - 1) * beta%y * beta%y / beta_sq
    L(3, 4) = (gamma - 1) * beta%y * beta%z / beta_sq

    L(4, 1) = L(1, 4)
    L(4, 2) = L(2, 4)
    L(4, 2) = L(3, 4)
    L(4, 4) = 1._8 + (gamma - 1) * beta%z * beta%z / beta_sq

    ! Prepare initial four-vector
    F1 = [1._8, d%x, d%y, d%z] * nu

    ! Compute final four-vector
    F2 = matmul(L, F1)

    ! Extract required quantities
    nu = F2(1)
    d%x = F2(2) / nu
    d%y = F2(3) / nu
    d%z = F2(4) / nu

    ! Normalize output direction vector
    norm = sqrt(d%x * d%x + d%y * d%y + d%z * d%z)
    d%x = d%x / norm
    d%y = d%y / norm
    d%z = d%z / norm

    ! TODO - optimize this because we are wasting time converting a to v and vice-versa
    call vector3d_to_angle3d(d, a)

  end subroutine lorentz_boost

end module lorentz
