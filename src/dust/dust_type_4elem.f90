module type_dust

  use core_lib
  use mpi_hdf5_io

  implicit none
  save

  private
  public :: dust
  public :: dust_setup
  public :: dust_emit
  public :: dust_emit_peeloff
  public :: dust_scatter
  public :: dust_scatter_peeloff
  public :: dust_sample_emit_probability
  public :: dust_sample_j_nu
  public :: dust_sample_b_nu
  public :: dust_jnu_var_pos_frac

  public :: get_j_nu_interp
  public :: get_j_nu_binned
  public :: get_chi_nu_interp
  public :: get_chi_nu_binned

  type dust

     ! Metadata
     integer :: version

     ! Sublimation
     integer :: sublimation_mode = 0 ! 0=no, 1=fast, 2=slow, 3=capped
     real(dp) :: sublimation_specific_energy

     ! Optical properties
     integer :: n_nu                              ! number of frequencies
     real(dp),allocatable :: nu(:),log_nu(:)      ! Frequency
     real(dp),allocatable :: albedo_nu(:)         ! Albedo
     real(dp),allocatable :: chi_nu(:)          ! Opacity
     real(dp),allocatable :: kappa_nu(:)      ! Opacity to absorption

     ! Scattering matrix
     integer :: n_mu                              ! number of cos(theta)
     real(dp),allocatable :: mu(:)                ! Values of cos(theta) that P is tabulated for
     real(dp),allocatable :: P1(:,:),P2(:,:),P3(:,:),P4(:,:) ! 4-element P matrix
     real(dp),allocatable :: P1_cdf(:,:),P2_cdf(:,:),P3_cdf(:,:),P4_cdf(:,:) ! 4-element P matrix
     real(dp),allocatable :: I_max(:)
     real(dp) :: mu_min, mu_max                   ! range of mu values

     ! Mean opacities
     integer :: n_e                              ! Number of energies
     real(dp),allocatable :: specific_energy(:) ! Energy absorbed per unit mass
     real(dp),allocatable :: log10_specific_energy(:) ! Energy absorbed per unit mass [Log10]
     real(dp),allocatable :: chi_planck(:)       ! Planck mean opacity
     real(dp),allocatable :: kappa_planck(:)     ! Planck mean absoptive opacity
     real(dp),allocatable :: chi_inv_planck(:)    ! Rosseland mean opacity
     real(dp),allocatable :: kappa_inv_planck(:)  ! Rosseland mean opacity
     real(dp),allocatable :: chi_rosseland(:)    ! Rosseland mean opacity
     real(dp),allocatable :: kappa_rosseland(:)  ! Rosseland mean opacity

     ! Emissivity
     integer :: n_jnu                             ! number of emissivities
     character(len=1) :: emiss_var                ! type of independent emissivity variable
     real(dp),allocatable :: j_nu_var(:)          ! independent emissivity variable
     real(dp),allocatable :: log10_j_nu_var(:)    ! independent emissivity variable [Log10]
     type(pdf_dp),allocatable :: j_nu(:)          ! emissivity
     type(pdf_dp),allocatable :: b_nu(:)          ! emissivity divided by opacity (blackbodies for LTE dust)

     logical :: is_lte ! Whether the emissivities assume therma emission from LTE dust

     logical :: zero_p2  ! Whether the P2 term is zero (since it makes the mu sampling simpler)

  end type dust

contains

  subroutine dust_setup(group,d,beta)

    implicit none

    integer(hid_t),intent(in) :: group
    type(dust),intent(out)    :: d
    real(dp),intent(in)       :: beta
    integer :: i,j
    real(dp),allocatable :: emiss_nu(:), emiss_jnu(:,:)
    real(dp) :: norm, dmu
    character(len=100) :: path
    character(len=4) :: sublimation
    type(version) :: python_version
    integer :: dust_version

    ! Read dust file

    if(mp_exists_keyword(group, '.', 'python_version')) then
       call mp_read_keyword(group, '.', 'python_version', python_version%string)
       if(python_version < version('0.8.7')) then
          call error("setup_initial", "cannot read dust files made with the Python module before version 0.8.7")
       end if
    else
       call error("setup_initial", "cannot read dust files made with the Python module before version 0.8.7")
    end if

    call mp_read_keyword(group, '.', 'version', d%version)

    call mp_read_keyword(group, '.', 'emissvar', d%emiss_var)

    call mp_read_keyword(group, '.', 'lte', d%is_lte)

    if(d%emiss_var /= 'E') stop "Only emissvar='E' supported at this time"

    ! DUST SUBLIMATION

    call mp_read_keyword(group, '.', 'sublimation_mode', sublimation)

    select case(trim(sublimation))
    case('no')
       d%sublimation_mode = 0
    case('fast')
       d%sublimation_mode = 1
       call mp_read_keyword(group, '.', 'sublimation_specific_energy', d%sublimation_specific_energy)
    case('slow')
       d%sublimation_mode = 2
       call mp_read_keyword(group, '.', 'sublimation_specific_energy', d%sublimation_specific_energy)
    case('cap')
       d%sublimation_mode = 3
       call mp_read_keyword(group, '.', 'sublimation_specific_energy', d%sublimation_specific_energy)
    case default
       call error('setup_initial','Unknown dust sublimation mode: '//trim(sublimation))
    end select

    ! OPTICAL PROPERTIES

    path = 'optical_properties'
    call mp_table_read_column_auto(group,path,'nu',d%nu)
    call mp_table_read_column_auto(group,path,'albedo',d%albedo_nu)
    call mp_table_read_column_auto(group,path,'chi',d%chi_nu)
    call mp_table_read_column_auto(group,path,'P1',d%P1)
    call mp_table_read_column_auto(group,path,'P2',d%P2)
    call mp_table_read_column_auto(group,path,'P3',d%P3)
    call mp_table_read_column_auto(group,path,'P4',d%P4)

    ! Check whether there is an (I,Q) cross-term
    d%zero_p2 = all(d%P2 == 0._dp)

    ! Check for NaN values
    if(any(is_nan(d%nu))) call error("dust_setup","nu array contains NaN values")
    if(any(is_nan(d%albedo_nu))) call error("dust_setup","albedo_nu array contains NaN values")
    if(any(is_nan(d%chi_nu))) call error("dust_setup","chi_nu array contains NaN values")
    if(any(is_nan(d%P1))) call error("dust_setup","P1 matrix contains NaN values")
    if(any(is_nan(d%P2))) call error("dust_setup","P2 matrix contains NaN values")
    if(any(is_nan(d%P3))) call error("dust_setup","P3 matrix contains NaN values")
    if(any(is_nan(d%P4))) call error("dust_setup","P4 matrix contains NaN values")

    ! Find number of frequencies
    d%n_nu = size(d%nu)

    ! Compute log[nu]
    allocate(d%log_nu(d%n_nu))
    d%log_nu = log10(d%nu)

    ! Compute opacity to absorption
    allocate(d%kappa_nu(d%n_nu))
    d%kappa_nu = d%chi_nu * (1._dp - d%albedo_nu)

    ! Compute maximum scattering intensity vs wavelength
    allocate(d%I_max(d%n_nu))
    do j=1,d%n_nu
       d%I_max(j) = maxval(d%P1(:,j)+abs(d%P2(:,j)))
    end do

    path = 'scattering_angles'
    call mp_table_read_column_auto(group,path,'mu',d%mu)

    ! Check for NaN values
    if(any(is_nan(d%mu))) call error("dust_setup","mu array contains NaN values")

    ! Find number of scattering angles
    d%n_mu = size(d%mu)

    ! Find min and max
    d%mu_min = d%mu(1)
    d%mu_max = d%mu(d%n_mu)

    dmu = d%mu_max - d%mu_min

    ! Normalize scattering matrix. The probability distribution functions
    ! are normalized so that their total integrals are 4*pi (not 1)
    do j=1,d%n_nu
       norm = integral_linlog(d%mu, d%P1(:,j))
       if(norm.eq.0._dp) call error("dust_setup", "P1 matrix normalization is zero")
       d%P1(:,j) = d%P1(:,j) / norm * dmu
       d%P2(:,j) = d%P2(:,j) / norm * dmu
       d%P3(:,j) = d%P3(:,j) / norm * dmu
       d%P4(:,j) = d%P4(:,j) / norm * dmu
    end do

    ! Allocate cumulative scattering matrix elements
    allocate(d%P1_cdf(size(d%P1,1), size(d%P1,2)))
    allocate(d%P2_cdf(size(d%P2,1), size(d%P2,2)))
    allocate(d%P3_cdf(size(d%P3,1), size(d%P3,2)))
    allocate(d%P4_cdf(size(d%P4,1), size(d%P4,2)))

    ! Find cumulative scattering matrix elements
    do j=1,d%n_nu
       d%P1_cdf(:,j) = cumulative_integral(d%mu, d%P1(:,j))
       d%P2_cdf(:,j) = cumulative_integral(d%mu, d%P2(:,j))
       d%P3_cdf(:,j) = cumulative_integral(d%mu, d%P3(:,j))
       d%P4_cdf(:,j) = cumulative_integral(d%mu, d%P4(:,j))
       if(.not.all(d%P1_cdf(:,j)==0.)) d%P1_cdf(:,j) = d%P1_cdf(:,j) / d%P1_cdf(d%n_mu, j)
       if(.not.all(d%P2_cdf(:,j)==0.)) d%P2_cdf(:,j) = d%P2_cdf(:,j) / d%P2_cdf(d%n_mu, j)
       if(.not.all(d%P3_cdf(:,j)==0.)) d%P3_cdf(:,j) = d%P3_cdf(:,j) / d%P3_cdf(d%n_mu, j)
       if(.not.all(d%P4_cdf(:,j)==0.)) d%P4_cdf(:,j) = d%P4_cdf(:,j) / d%P4_cdf(d%n_mu, j)
    end do

    ! MEAN OPACITIES

    path = 'mean_opacities'
    call mp_table_read_column_auto(group,path,'specific_energy',d%specific_energy)
    call mp_table_read_column_auto(group,path,'chi_planck',d%chi_planck)
    call mp_table_read_column_auto(group,path,'kappa_planck',d%kappa_planck)
    call mp_table_read_column_auto(group,path,'chi_rosseland',d%chi_rosseland)
    call mp_table_read_column_auto(group,path,'kappa_rosseland',d%kappa_rosseland)

    ! In version 1 there was a bug that caused the Rosseland mean opacity to
    ! be mis-computed (it was in fact computing the Planck inverse opacity).
    ! However, this only affects models that use the PDA, so we only need to
    ! raise an error for these. For the MRW, what was needed was *actually*
    ! the reciprocal Planck opacity so models using this are fine. So we
    ! don't raise an error here but we keep track of the version in the dust
    ! object and can raise an error in the setup part of the code.

    if(d%version == 1) then
       call mp_table_read_column_auto(group,path,'chi_rosseland',d%chi_inv_planck)
       call mp_table_read_column_auto(group,path,'kappa_rosseland',d%kappa_inv_planck)
    else
       call mp_table_read_column_auto(group,path,'chi_inv_planck',d%chi_inv_planck)
       call mp_table_read_column_auto(group,path,'kappa_inv_planck',d%kappa_inv_planck)
    end if

    ! Check for NaN values
    if(any(is_nan(d%specific_energy))) call error("dust_setup","specific_energy array contains NaN values")
    if(any(is_nan(d%chi_planck))) call error("dust_setup","chi_planck array contains NaN values")
    if(any(is_nan(d%kappa_planck))) call error("dust_setup","kappa_planck array contains NaN values")
    if(any(is_nan(d%chi_inv_planck))) call error("dust_setup","chi_inv_planck array contains NaN values")
    if(any(is_nan(d%kappa_inv_planck))) call error("dust_setup","kappa_inv_planck array contains NaN values")
    if(any(is_nan(d%chi_rosseland))) call error("dust_setup","chi_planck array contains NaN values")
    if(any(is_nan(d%kappa_rosseland))) call error("dust_setup","kappa_rosseland array contains NaN values")

    d%n_e = size(d%specific_energy)
    allocate(d%log10_specific_energy(d%n_e))
    d%log10_specific_energy = log10(d%specific_energy)

    ! Check that specific energy is monotically increasing (important for interpolation)
    do i=2,d%n_e
       if(d%specific_energy(i) < d%specific_energy(i-1)) then
          call error("dust_setup","energy per unit mass is not monotonically increasing")
       end if
    end do

    ! need to check monotonically increases

    ! EMISSIVITIES

    path = 'emissivities'
    call mp_table_read_column_auto(group,path,'nu',emiss_nu)
    call mp_table_read_column_auto(group,path,'jnu',emiss_jnu)

    ! Check for NaN values
    if(any(is_nan(emiss_nu))) call error("dust_setup","emiss_nu array contains NaN values")
    if(any(is_nan(emiss_jnu))) call error("dust_setup","emiss_jnu array contains NaN values")

    path = 'emissivity_variable'
    select case(d%emiss_var)
    case('E')
       call mp_table_read_column_auto(group,path,'specific_energy',d%j_nu_var)
       if(any(is_nan(d%j_nu_var))) call error("dust_setup","emissivity variable array contains NaN values")
    end select

    ! Find number of emissivites
    d%n_jnu = size(d%j_nu_var)
    allocate(d%log10_j_nu_var(d%n_jnu))
    d%log10_j_nu_var = log10(d%j_nu_var)

    ! Allocate emissivity PDF
    allocate(d%j_nu(d%n_jnu))
    allocate(d%b_nu(d%n_jnu))

    do i=1,d%n_jnu
       call set_pdf(d%j_nu(i),emiss_nu,emiss_jnu(i,:),log=.true.)
       call set_pdf(d%b_nu(i),emiss_nu,emiss_jnu(i,:) / interp1d_loglog(d%nu, d%kappa_nu, emiss_nu),log=.true.)
    end do

  end subroutine dust_setup

  subroutine dust_jnu_var_pos_frac(d,specific_energy,jnu_var_id,jnu_var_frac)
    implicit none
    type(dust),intent(in) :: d
    real(dp),intent(in) :: specific_energy
    integer,intent(out) :: jnu_var_id
    real(dp),intent(out) :: jnu_var_frac
    real(dp) :: jnu_var

    select case(d%emiss_var)
    case('E')
       jnu_var = specific_energy
    end select

    if(jnu_var < d%j_nu_var(1)) then
       jnu_var_id = 1
       jnu_var_frac = 0._dp
    else if(jnu_var > d%j_nu_var(size(d%j_nu_var))) then
       jnu_var_id = size(d%j_nu_var) - 1
       jnu_var_frac = 1._dp
    else
       jnu_var_id = locate(d%j_nu_var,jnu_var)
       jnu_var_frac = (log10(jnu_var) - d%log10_j_nu_var(jnu_var_id)) &
            &       / (d%log10_j_nu_var(jnu_var_id + 1) - d%log10_j_nu_var(jnu_var_id))
    end if

  end subroutine dust_jnu_var_pos_frac

  subroutine dust_emit_peeloff(d,nu,a,s,a_req)
    implicit none
    type(dust),intent(in)          :: d
    real(dp),intent(in)            :: nu
    type(angle3d_dp),intent(inout) :: a
    type(stokes_dp),intent(inout)  :: s
    type(angle3d_dp),intent(in)    :: a_req
    ! The probability distribution function for the redistribution is
    ! normalized so that its total integral is 4*pi (not 1)
    a = a_req
  end subroutine dust_emit_peeloff

  subroutine dust_emit(d,jnu_var_id,jnu_var_frac,nu,a,s,energy_scaling)

    implicit none

    type(dust),intent(in)          :: d
    integer,intent(in)             :: jnu_var_id
    real(dp),intent(in)            :: jnu_var_frac
    type(angle3d_dp),intent(out)   :: a
    type(stokes_dp),intent(out)    :: s
    real(dp),intent(out)           :: nu
    real(dp),intent(out)           :: energy_scaling

    call dust_sample_j_nu(d,jnu_var_id,jnu_var_frac,nu)

    s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    energy_scaling = 1.

    call random_sphere_angle3d(a)

  end subroutine dust_emit

  subroutine dust_sample_emit_probability(d,jnu_var_id,jnu_var_frac,nu, prob)

    implicit none

    type(dust),intent(in)          :: d
    integer,intent(in)             :: jnu_var_id
    real(dp),intent(in)            :: jnu_var_frac, nu
    real(dp),intent(out)           :: prob

    real(dp) :: prob1,prob2

    prob1 = interpolate_pdf(d%j_nu(jnu_var_id), nu, bounds_error=.false., fill_value=0._dp)
    prob2 = interpolate_pdf(d%j_nu(jnu_var_id+1), nu, bounds_error=.false., fill_value=0._dp)

    if(prob1.eq.0._dp.or.prob2.eq.0._dp) then
       prob = 0._dp
    else
       prob = log10(prob1) + jnu_var_frac * (log10(prob2) - log10(prob1))
       prob = 10._dp**prob
    end if

  end subroutine dust_sample_emit_probability

  subroutine dust_sample_j_nu(d,jnu_var_id,jnu_var_frac,nu)

    implicit none

    type(dust),intent(in)          :: d
    integer,intent(in)             :: jnu_var_id
    real(dp),intent(in)            :: jnu_var_frac
    real(dp),intent(out)           :: nu

    real(dp) :: nu1,nu2,xi

    call random(xi)

    nu1 = sample_pdf(d%j_nu(jnu_var_id),xi)
    nu2 = sample_pdf(d%j_nu(jnu_var_id+1),xi)

    nu = log10(nu1) + jnu_var_frac * (log10(nu2) - log10(nu1))
    nu = 10._dp**nu

  end subroutine dust_sample_j_nu

  subroutine dust_sample_b_nu(d,jnu_var_id,jnu_var_frac,nu)

    implicit none

    type(dust),intent(in)          :: d
    integer,intent(in)             :: jnu_var_id
    real(dp),intent(in)            :: jnu_var_frac
    real(dp),intent(out)           :: nu

    real(dp) :: nu1,nu2,xi

    call random(xi)

    nu1 = sample_pdf(d%b_nu(jnu_var_id),xi)
    nu2 = sample_pdf(d%b_nu(jnu_var_id+1),xi)

    nu = log10(nu1) + jnu_var_frac * (log10(nu2) - log10(nu1))
    nu = 10._dp**nu

  end subroutine dust_sample_b_nu

  subroutine dust_scatter_peeloff(d,nu,a,s,a_req)
    implicit none
    type(dust),intent(in)          :: d
    real(dp),intent(in)            :: nu
    type(angle3d_dp),intent(inout) :: a
    type(stokes_dp),intent(inout)  :: s
    type(angle3d_dp),intent(in)    :: a_req
    type(angle3d_dp) :: a_scat
    real(dp) :: P1,P2,P3,P4
    call difference_angle3d(a, a_req, a_scat)
    if(a_scat%cost < d%mu_min .or. a_scat%cost > d%mu_max) then
       s%i = 0.
       s%q = 0.
       s%u = 0.
       s%v = 0.
    else
       P1 = interp2d(d%mu,d%nu,d%P1,a_scat%cost,nu)
       P2 = interp2d(d%mu,d%nu,d%P2,a_scat%cost,nu)
       P3 = interp2d(d%mu,d%nu,d%P3,a_scat%cost,nu)
       P4 = interp2d(d%mu,d%nu,d%P4,a_scat%cost,nu)
       call scatter_stokes(s,a,a_scat,a_req,P1,P2,P3,P4)
    end if
    a = a_req
  end subroutine dust_scatter_peeloff

  subroutine dust_scatter(d,nu,a,s)

    implicit none

    type(dust),intent(in)                :: d
    type(angle3d_dp),intent(inout)       :: a
    type(stokes_dp),intent(inout)        :: s

    real(dp) :: nu

    type(angle3d_dp) :: a_scat
    type(angle3d_dp) :: a_final

    real(dp) :: P1,P2,P3,P4,norm

    real(dp) :: c1, c2, ctot, cdf1, cdf2, xi
    real(dp) :: sin_2_i1,cos_2_i1

    integer :: imin, imax, imu, inu

    integer :: iter
    integer,parameter :: maxiter = 1000000

    !#############################################################################
    !
    ! In order to sample the scattering angle, we first sample two angles
    ! theta and phi uniformly.
    !
    ! We then calculate the new value of I using these values, and the previous
    ! values of the Stokes parameters, and we decide whether to keep it using
    ! the rejection criterion
    !
    !#############################################################################

    call random_sphere_angle3d(a_scat)

    sin_2_i1 =         2._dp * a_scat%sinp * a_scat%cosp
    cos_2_i1 = 1._dp - 2._dp * a_scat%sinp * a_scat%sinp

    c1 = s%I
    c2 = (cos_2_i1 * s%Q - sin_2_i1 * s%U)
    ctot = c1 + c2
    c1 = c1 / ctot
    c2 = c2 / ctot

    imin = 1
    imax = d%n_mu
    inu = locate(d%nu, nu)
    ! TODO: interpolate in nu as well

    if(inu==-1) then

       ! Frequency is out of bounds, use isotropic scattering
       P1 = 1._dp
       P2 = 0._dp
       P3 = 1._dp
       P4 = 0._dp

    else

       call random(xi)

       if(d%zero_p2) then  ! no I and Q cross term, simple sampling
          do iter=1,maxiter
             imu = (imax + imin) / 2
             cdf1 = d%P1_cdf(imu, inu)
             cdf2 = d%P1_cdf(imu+1, inu)
             if(xi > cdf2) then
                imin = imu
             else if(xi < cdf1) then
                imax = imu
             else
                exit
             end if
             if(imin==imax) stop "ERROR: in sampling mu for scattering"
          end do
       else
          do iter=1,maxiter
             imu = (imax + imin) / 2
             cdf1 = c1 * d%P1_cdf(imu, inu) + c2 * d%P2_cdf(imu, inu)
             cdf2 = c1 * d%P1_cdf(imu+1, inu) + c2 * d%P2_cdf(imu+1, inu)
             if(xi > cdf2) then
                imin = imu
             else if(xi < cdf1) then
                imax = imu
             else
                exit
             end if
             if(imin==imax) stop "ERROR: in sampling mu for scattering"
          end do
       end if

       if(iter==maxiter+1) stop "ERROR: stuck in do loop in dust_scatter"

       a_scat%cost = (xi - cdf1) / (cdf2 - cdf1) * (d%mu(imu+1) - d%mu(imu)) + d%mu(imu)
       a_scat%sint = sqrt(1._dp - a_scat%cost*a_scat%cost)

       P1 = interp2d(d%mu,d%nu,d%P1,a_scat%cost,nu)
       P2 = interp2d(d%mu,d%nu,d%P2,a_scat%cost,nu)
       P3 = interp2d(d%mu,d%nu,d%P3,a_scat%cost,nu)
       P4 = interp2d(d%mu,d%nu,d%P4,a_scat%cost,nu)

    end if

    ! Find new photon direction
    call rotate_angle3d(a_scat,a,a_final)

    ! Compute how the stokes parameters are changed by the interaction
    call scatter_stokes(s,a,a_scat,a_final,P1,P2,P3,P4)

    ! Change photon direction
    a = a_final

    norm = 1._dp / S%I

    S%I = 1._dp
    S%Q = S%Q * norm
    S%U = S%U * norm
    S%V = S%V * norm

  end subroutine dust_scatter

  !#############################################################################
  !
  ! To find how the stokes parameters S = (I,Q,U,V) change with the scattering
  ! interaction, use the following equation:
  !
  ! S = L( pi - i_2 ) * R * L( - i_1 ) * S'
  !
  ! S' is the old set of Stokes parameters
  ! L ( - i_1 ) is a rotation matrix to rotate into the plane of scattering
  ! R calculates the scattering function
  ! L ( pi - i_2 ) rotates back to the observer's frame of reference
  !
  ! The rotation matrix L is given by
  !
  !          /  1  |      0      |      0      |  0  \
  ! L(psi) = |  0  | +cos(2*psi) | +sin(2*psi) |  0  |
  !          |  0  | -sin(2*psi) | +cos(2*psi) |  0  |
  !          \  0  |      0      |      0      |  1  /
  !
  ! The scattering matrix can have various number of elements.
  !
  ! The electron or dust scattering properties are recorded in the R matrix.
  !
  ! For example, a four element matrix could be:
  !
  !                /  P1  P2  0   0  \
  ! R(theta) = a * |  P2  P1  0   0  |
  !                |  0   0   P3 -P4 |
  !                \  0   0   P4  P3 /
  !
  ! The values of P1->4 can either be found from an analytical function, or
  ! read in from files.
  !
  !#############################################################################

  subroutine scatter_stokes(s,a_coord,a_scat,a_final,P1,P2,P3,P4)

    implicit none

    type(angle3d_dp),intent(in)    :: a_coord     ! The photon direction angle
    type(angle3d_dp),intent(in)    :: a_scat      ! The photon scattering angle
    type(angle3d_dp),intent(in)    :: a_final     ! The final photon direction
    type(stokes_dp),intent(inout)  :: s           ! The Stokes parameters of the photon
    real(dp),intent(in)            :: P1,P2,P3,P4 ! 4-element matrix elements

    ! Spherical trigonometry
    real(dp) :: cos_a,sin_a
    real(dp) :: cos_b,sin_b
    real(dp) :: cos_c,sin_c
    real(dp) :: cos_big_a,sin_big_a
    real(dp) :: cos_big_b,sin_big_b
    real(dp) :: cos_big_c,sin_big_c

    ! Local
    real(dp) :: cos_i2,cos_2_i2
    real(dp) :: sin_i2,sin_2_i2
    real(dp) :: cos_2_alpha,cos_2_beta
    real(dp) :: sin_2_alpha,sin_2_beta
    real(dp) :: RLS1,RLS2,RLS3,RLS4

    ! The general spherical trigonometry routines in type_angle3d have served
    ! us well this far, but now we need to compute a specific angle in the
    ! spherical triangle. The meaning of the angles is as follows:

    ! a =   old theta angle (initial direction angle)
    ! b = local theta angle (scattering or emission angle)
    ! c =   new theta angle (final direction angle)

    ! A = what we want to calculate here
    ! B = new phi - old phi
    ! C = local phi angle (scattering or emission angle)

    cos_a = a_coord%cost
    sin_a = a_coord%sint

    cos_b = a_scat%cost
    sin_b = a_scat%sint

    cos_c = a_final%cost
    sin_c = a_final%sint

    cos_big_b = a_coord%cosp * a_final%cosp + a_coord%sinp * a_final%sinp
    sin_big_b = a_coord%sinp * a_final%cosp - a_coord%cosp * a_final%sinp

    cos_big_C = a_scat%cosp
    sin_big_C = abs(a_scat%sinp)

    if(sin_big_c < 10. * tiny(1._dp) .and. sin_c < 10. * tiny(1._dp)) then
       cos_big_a = - cos_big_b * cos_big_c
       sin_big_a = sqrt(1._8 - cos_big_a * cos_big_a)
    else
       cos_big_a = (cos_a - cos_b * cos_c) / (sin_b * sin_c)
       sin_big_a = + sin_big_c * sin_a / sin_c
    end if

    cos_i2 = cos_big_a
    sin_i2 = sin_big_a

    cos_2_i2    = 1._dp - 2._dp * sin_i2 * sin_i2
    sin_2_i2    =         2._dp * sin_i2 * cos_i2

    cos_2_alpha = 1._dp - 2._dp * a_scat%sinp * a_scat%sinp
    sin_2_alpha =       - 2._dp * a_scat%sinp * a_scat%cosp

    if(a_scat%sinp < 0.) then
       cos_2_beta =  cos_2_i2
       sin_2_beta =  sin_2_i2
    else
       cos_2_beta =  cos_2_i2
       sin_2_beta = -sin_2_i2
    end if

    RLS1 =   P1 * S%I + P2 * ( + cos_2_alpha * S%Q + sin_2_alpha * S%U )
    RLS2 =   P2 * S%I + P1 * ( + cos_2_alpha * S%Q + sin_2_alpha * S%U )
    RLS3 = - P4 * S%V + P3 * ( - sin_2_alpha * S%Q + cos_2_alpha * S%U )
    RLS4 =   P3 * S%V + P4 * ( - sin_2_alpha * S%Q + cos_2_alpha * S%U )

    S%I = RLS1
    S%Q = + cos_2_beta * RLS2 + sin_2_beta * RLS3
    S%U = - sin_2_beta * RLS2 + cos_2_beta * RLS3
    S%V = RLS4

  end subroutine scatter_stokes

  elemental real(dp) function B_nu(nu,T)
    implicit none
    real(dp),intent(in) :: nu,T
    real(dp),parameter :: a = two * h_cgs / c_cgs / c_cgs
    real(dp),parameter :: b = h_cgs / k_cgs
    B_nu = a * nu * nu * nu / ( exp(b*nu/T) - one)
  end function B_nu

  elemental real(dp) function dB_nu_over_dT(nu,T)
    implicit none
    real(dp),intent(in) :: nu,T
    real(dp),parameter :: a = two * h_cgs * h_cgs / c_cgs / c_cgs / k_cgs
    real(dp),parameter :: b = h_cgs / k_cgs
    dB_nu_over_dT = a * nu * nu * nu * nu * exp(b*nu/T) / (exp(b*nu/T) - one)**2.
  end function dB_nu_over_dT

  function get_j_nu_interp(d, nu, jnu_var_id) result(j_nu)

    implicit none

    type(dust),intent(in) :: d
    real(dp),intent(in) :: nu(:)
    integer,intent(in) :: jnu_var_id
    real(dp) :: j_nu(size(nu))

    j_nu = interp1d_loglog(d%j_nu(jnu_var_id)%x, d%j_nu(jnu_var_id)%pdf, &
         &                 nu, bounds_error=.false., fill_value=0._dp)

  end function get_j_nu_interp

  function get_j_nu_binned(d, n_nu, nu_min, nu_max, jnu_var_id) result(j_nu)

    implicit none

    type(dust),intent(in) :: d
    integer,intent(in) :: n_nu
    real(dp),intent(in) :: nu_min, nu_max
    integer,intent(in) :: jnu_var_id
    real(dp) :: j_nu(n_nu)

    integer :: inu
    real(dp) :: numin, numax
    real(dp) :: log10_nu_min, log10_nu_max

    log10_nu_min = log10(nu_min)
    log10_nu_max = log10(nu_max)

    do inu=1, n_nu

       numin = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu - 1, dp) / real(n_nu, dp))
       numax = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu, dp) / real(n_nu, dp))

       j_nu(inu) = integral_loglog(d%j_nu(jnu_var_id)%x, d%j_nu(jnu_var_id)%pdf, numin, numax)

    end do

    j_nu = j_nu / integral_loglog(d%j_nu(jnu_var_id)%x, d%j_nu(jnu_var_id)%pdf)

  end function get_j_nu_binned

  function get_chi_nu_interp(d, nu) result(chi_nu)

    implicit none

    type(dust),intent(in) :: d
    real(dp),intent(in) :: nu(:)
    real(dp) :: chi_nu(size(nu))

    chi_nu = interp1d_loglog(d%nu, d%chi_nu, &
         &                   nu, bounds_error=.false., fill_value=0._dp)

  end function get_chi_nu_interp

  function get_chi_nu_binned(d, n_nu, nu_min, nu_max) result(chi_nu)

    implicit none

    type(dust),intent(in) :: d
    integer,intent(in) :: n_nu
    real(dp),intent(in) :: nu_min, nu_max
    real(dp) :: chi_nu(n_nu)

    integer :: inu
    real(dp) :: numin, numax
    real(dp) :: log10_nu_min, log10_nu_max

    log10_nu_min = log10(nu_min)
    log10_nu_max = log10(nu_max)

    do inu=1, n_nu

       numin = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu - 1, dp) / real(n_nu, dp))
       numax = 10._dp**(log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu, dp) / real(n_nu, dp))

       chi_nu(inu) = integral_loglog(d%nu, d%chi_nu, numin, numax) / (numax - numin)

    end do

  end function get_chi_nu_binned

end module type_dust
