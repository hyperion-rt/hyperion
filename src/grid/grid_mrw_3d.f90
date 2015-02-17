module grid_mrw

  ! Optimization: could pre-compute diffusion coefficient just for masked (valid) cells

  use core_lib
  use type_photon
  use grid_geometry
  use grid_physics
  use dust_main, only : kappa_planck, chi_inv_planck, n_dust, d
  use type_dust, only : dust_sample_b_nu

  implicit none
  save

  private
  public :: prepare_mrw
  public :: grid_do_mrw
  public :: grid_do_mrw_noenergy
  public :: mrw_peeloff

  logical :: debug = .false.

  integer,parameter :: ncdf = 100
  real(dp) :: xcdf(ncdf), ycdf(ncdf)

  real(dp),allocatable :: diff_coeff(:)

contains

  subroutine prepare_mrw()

    implicit none

    integer :: ic, id

    real(dp) :: total_alpha_inv_planck
    ! eta_nu values have been read into dust type

    call initialize_cumulative()

    if(.not.allocated(diff_coeff)) allocate(diff_coeff(geo%n_cells))

    diff_coeff = 0._dp

    ! loop over cells
    do ic=1,geo%n_cells
       total_alpha_inv_planck = 0._dp
       do id=1,n_dust
          total_alpha_inv_planck = total_alpha_inv_planck + &
               &density(ic, id) * chi_inv_planck(id, specific_energy(ic, id))
       end do
       diff_coeff(ic) = 1._dp / 3._dp / total_alpha_inv_planck
    end do

  end subroutine prepare_mrw

  subroutine grid_do_mrw(p)

    implicit none

    type(photon),intent(inout) :: p

    real(dp) :: R0
    real(dp) :: e,y,ct
    type(vector3d_dp) :: dr
    integer :: id

    ! Find distance to closest wall
    R0 = distance_to_closest_wall(p)

    ! Sample P0
    y = sample_cumulative()

    ! Solve (8) for ct
    ct = -log(y) / diff_coeff(p%icell%ic) * (R0/pi)**2.

    ! Deposit energy in cell
    do id=1,n_dust
       if(density(p%icell%ic, id) > 0._dp) then
          ! Insert ct into (9), get energy deposited for Lucy method
          e = p%energy * ct * kappa_planck(id, specific_energy(p%icell%ic, id))
          specific_energy_sum(p%icell%ic, id) = specific_energy_sum(p%icell%ic, id) + e
       end if
    end do

    ! Place photon on sphere surface, optionally sample new direction
    call random_sphere_vector3d(dr)

    p%r = p%r + dr * R0

    ! Save last direction before scattering/absorbing
    p%a_prev = p%a
    p%v_prev = p%v
    p%s_prev = p%s

    call random_sphere_angle3d(p%a)
    call angle3d_to_vector3d(p%a, p%v)

    id = select_dust_chi_rho(p)
    call dust_sample_b_nu(d(id), jnu_var_id(p%icell%ic, id), jnu_var_frac(p%icell%ic, id), p%nu)

    ! For peeloff, we are going to assume that the radiation is isotropic.
    ! This is not quite exact, but is not likely to matter much.
    p%last_isotropic = .true.
    p%dust_id = id
    p%last = 'de'

  end subroutine grid_do_mrw

  subroutine grid_do_mrw_noenergy(p)

    implicit none

    type(photon),intent(inout) :: p

    real(dp) :: R0
    type(vector3d_dp) :: dr
    integer :: id

    ! Find distance to closest wall
    R0 = distance_to_closest_wall(p)

    ! Place photon on sphere surface, optionally sample new direction
    call random_sphere_vector3d(dr)

    p%r = p%r + dr * R0

    ! Save last direction before scattering/absorbing
    p%a_prev = p%a
    p%v_prev = p%v
    p%s_prev = p%s

    ! Sample frequency from eta_nu
    call random_sphere_angle3d(p%a)
    call angle3d_to_vector3d(p%a, p%v)

    id = select_dust_chi_rho(p)
    call dust_sample_b_nu(d(id), jnu_var_id(p%icell%ic, id), jnu_var_frac(p%icell%ic, id), p%nu)

    ! For peeloff, we are going to assume that the radiation is isotropic.
    ! This is not quite exact, but is not likely to matter much.
    p%last_isotropic = .true.
    p%dust_id = id
    p%last = 'me'

  end subroutine grid_do_mrw_noenergy

  subroutine mrw_peeloff(p,a_req)
    implicit none
    type(photon),intent(inout) :: p
    type(angle3d_dp),intent(in)    :: a_req
    p%a = a_req
  end subroutine mrw_peeloff

  subroutine initialize_cumulative()

    ! The purpose of this function is to compute equation (6) of Min et al (2009):
    !
    ! P(t) = 2 * Sum_{n=1}^{infinity} (-1)^(n+1) * y^(n^2)
    !
    ! This function requires high n values close to y=1, so it is best to pre-compute it and to then interpolate

    implicit none

    integer(idp) :: i,j
    real(dp) :: term

    ycdf = 0._dp

    do i=1,ncdf

       xcdf(i) = real(i-1, dp)/real(ncdf-1, dp)

       if(i==ncdf) then
          ycdf(i) = 0.5_dp
       else
          j = 0
          do
             j = j + 1
             term = xcdf(i)**(j**2)
             if(term == 0._dp) exit
             if(mod(j, 2_idp)==0) then
                ycdf(i) = ycdf(i) - term
             else
                ycdf(i) = ycdf(i) + term
             end if
          end do
       end if
    end do

    ycdf = ycdf * 2._dp

  end subroutine initialize_cumulative

  real(dp) function sample_cumulative() result(xi)
    implicit none
    call random(xi)
    xi = interp1d(ycdf, xcdf, xi)
  end function sample_cumulative

end module grid_mrw
