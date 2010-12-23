module grid_physics

  use core_lib
  use type_photon
  use type_grid_cell
  use grid_io, only : read_grid_4d, grid_exists
  use dust_main ! many variables and routines
  use type_dust
  use grid_geometry
  use settings
  use mpi_core

  implicit none
  save

  private
  public :: setup_grid_physics
  public :: update_temperature
  public :: sublimate_dust
  public :: update_alpha_rosseland
  public :: update_dtau_rosseland
  public :: update_mean_temperature
  public :: update_energy_abs
  public :: update_energy_abs_tot
  public :: select_dust_chi_rho
  public :: select_dust_temp4_rho
  public :: emit_from_grid
  public :: precompute_jnu_var
  public :: tau_rosseland_to_closest_wall
  public :: specific_energy_abs_converged

  ! Density (immutable)
  real(dp),allocatable,target, public :: density(:,:)
  real(dp),allocatable, public :: density_original(:,:)

  ! Variable quantities (made public for MPI)
  real(dp),allocatable, public :: temperature(:,:)
  integer(idp),allocatable, public :: n_photons(:)
  integer(idp),allocatable, public :: last_photon_id(:)
  real(dp),allocatable, public :: specific_energy_abs(:,:)
  real(dp),allocatable, public :: energy_abs_tot(:)

  real(dp), allocatable,target, public :: temperature_mean(:)
  real(dp), allocatable,target, public :: dtau_rosseland(:,:)
  real(dp), allocatable,target, public :: alpha_rosseland(:)

  integer, allocatable, public :: jnu_var_id(:,:)
  real(dp), allocatable, public :: jnu_var_frac(:,:)

  logical, allocatable, public :: empty_neighbor(:)

  ! Temporary variables (used for convenience in external modules)
  real(dp), allocatable, public :: tmp_column_density(:)

  integer :: id

  logical :: debug = .false.

  type(pdf_discrete_dp) :: absorption

contains

  real(dp) function tau_rosseland_to_closest_wall(p) result(tau)
    implicit none
    type(photon),intent(in) :: p
    tau = alpha_rosseland(p%icell%ic) * distance_to_closest_wall(p)
  end function tau_rosseland_to_closest_wall

  integer function select_dust_chi_rho(p) result(id_select)
    implicit none
    type(photon),intent(in) :: p
    if(n_dust == 1) then
       id_select = 1
    else
       do id=1,n_dust
          absorption%pdf(id) = p%current_chi(id) * density(p%icell%ic, id)
       end do
       call find_cdf(absorption)
       id_select = sample_pdf(absorption)    
    end if
  end function select_dust_chi_rho

  integer function select_dust_temp4_rho(icell) result(id_select)
    implicit none
    type(grid_cell),intent(in) :: icell
    do id=1,n_dust
       if(is_lte_dust(id)) then ! temperature doesn't mean anything if not LTE dust
          absorption%pdf(id) = temperature(icell%ic, id)**4 * density(icell%ic, id)
       else
          absorption%pdf(id) = 0._dp
       end if
    end do
    call find_cdf(absorption)
    id_select = sample_pdf(absorption)    
  end function select_dust_temp4_rho

  subroutine setup_grid_physics(group, use_mrw, use_pda)

    implicit none

    integer(hid_t),intent(in) :: group
    logical,intent(in) :: use_mrw, use_pda

    ! Density
    allocate(density(geo%n_cells, n_dust))

    ! Temperature
    allocate(temperature(geo%n_cells, n_dust))

    if(n_dust > 0) then

       if(main_process()) write(*,'(" [grid_physics] reading density grid")')

       ! Read in density
       call read_grid_4d(group, 'Density', density, geo)

       ! Check number of dust types for density
       if(size(density, 2).ne.n_dust) call error("setup_grid","density array has wrong number of dust types")

       ! If density difference is requested, save original density
       if(output_density_diff.ne.'none') then
          allocate(density_original(geo%n_cells, n_dust))
          density_original = density
       end if

       if(grid_exists(group, 'Temperature')) then

          if(main_process()) write(*,'(" [grid_physics] reading temperature grid")')

          ! Read in temperature
          call read_grid_4d(group, 'Temperature', temperature, geo)

          ! Check number of dust types for temperature
          if(size(temperature, 2).ne.n_dust) call error("setup_grid","temperature array has wrong number of dust types")

          ! Check if any of the temperatures are below the minimum requested
          if(any(temperature < minimum_temperature)) then
             call warn("setup_grid_physics", "some of the initial temeperatures provided are below the requested minimum (resetting)")
             where(temperature < minimum_temperature)
                temperature = minimum_temperature
             end where
          end if

       else

          ! Set all temperatures to minimum requested
          temperature = minimum_temperature

       end if

    end if

    ! Column density for peeling-off
    allocate(tmp_column_density(n_dust))

    ! Specific energy absorbed
    allocate(specific_energy_abs(geo%n_cells, n_dust))
    specific_energy_abs = 0._dp

    ! Total energy absorbed
    allocate(energy_abs_tot(n_dust))
    energy_abs_tot = 0._dp

    ! Emissivity index and interpolation fraction
    allocate(jnu_var_id(geo%n_cells, n_dust))
    allocate(jnu_var_frac(geo%n_cells, n_dust))

    ! Modified Random Walk
    if(use_mrw) then

       ! Rosseland extinction coefficient
       allocate(alpha_rosseland(geo%n_cells))
       alpha_rosseland = 0._dp

    end if

    ! Partial Diffusion Approximation 
    if(use_pda) then

       ! Number of photons in each cell
       allocate(n_photons(geo%n_cells))
       n_photons = 0

       ! ID of last photon in each cell
       allocate(last_photon_id(geo%n_cells))
       last_photon_id = 0

       ! Rosseland optical depth across each cell (in 3 directions)
       allocate(dtau_rosseland(geo%n_cells, 3))
       dtau_rosseland = 0._dp

       ! Mean temperature across dust types
       allocate(temperature_mean(geo%n_cells))
       temperature_mean = minimum_temperature

    end if

    ! Create PDF for absorption in each cell
    call allocate_pdf(absorption,n_dust)

    ! Update energy absorbed in each cell to match temperature
    call update_energy_abs()

  end subroutine setup_grid_physics

  subroutine update_alpha_rosseland()

    implicit none

    integer :: ic

    if(main_process()) write(*,'(" [grid_physics] pre-computing Rosseland absorption coefficient")')
    alpha_rosseland = 0._dp

    do ic=1,geo%n_cells
       do id=1,n_dust
          alpha_rosseland(ic) = alpha_rosseland(ic) &
               & + density(ic,id) &
               & * chi_rosseland(id, temperature(ic,id))
       end do
    end do

  end subroutine update_alpha_rosseland

  subroutine update_dtau_rosseland()

    implicit none

    integer :: ic, idir

    if(main_process()) write(*,'(" [grid_physics] pre-computing Rosseland optical depth")')
    dtau_rosseland = 0._dp
    do ic=1,geo%n_cells
       do id=1,n_dust
          do idir=1,3
             dtau_rosseland(ic,idir) = dtau_rosseland(ic,idir) + &
                  &density(ic,id) * chi_rosseland(id, temperature(ic,id)) * geo%width(ic,idir)
          end do
       end do
    end do

  end subroutine update_dtau_rosseland

  subroutine update_mean_temperature()

    where(sum(density, dim=2) > 0.)
       temperature_mean = (sum(density * temperature**4., dim=2) / sum(density, dim=2))**0.25
    elsewhere
       temperature_mean = minimum_temperature
    end where

  end subroutine update_mean_temperature

  subroutine sublimate_dust()

    implicit none
    integer :: ic

    select case(dust_sublimation_mode)
    case(1)
       if(any(temperature > dust_sublimation_temperature)) then
          write(*,'(" [sublimate_dust] removing dust in ",I0," cells")') count(temperature > dust_sublimation_temperature)
          where(temperature > dust_sublimation_temperature)
             temperature = minimum_temperature
             density = 0.
          end where
       end if
    case(2)
       if(any(temperature > dust_sublimation_temperature)) then
          write(*,'(" [sublimate_dust] resetting density due to sublimation in ",I0," cells")') &
               & count(temperature > dust_sublimation_temperature)
          do ic=1,geo%n_cells
             do id=1,n_dust
                if (temperature(ic,id) > dust_sublimation_temperature) then
                   density(ic,id) = density(ic,id) &
                        & * (dust_sublimation_temperature / temperature(ic,id))**4 &
                        & * (kappa_planck(id, dust_sublimation_temperature) &
                        & / kappa_planck(id, temperature(ic,id))) &
                        & * (chi_rosseland(id, temperature(ic,id)) &
                        & / chi_rosseland(id, dust_sublimation_temperature))**2
                   temperature(ic,id) = dust_sublimation_temperature
                end if
             end do
          end do
       end if
    case(3)
       if(any(temperature > dust_sublimation_temperature)) then
          write(*,'(" [sublimate_dust] capping dust temperature in ",I0," cells")') count(temperature > dust_sublimation_temperature)
          where(temperature > dust_sublimation_temperature)
             temperature = dust_sublimation_temperature
          end where
       end if
    end select

    call update_energy_abs() ! Update energy absorbed in each cell to match temperature

  end subroutine sublimate_dust

  subroutine update_temperature()

    implicit none

    integer :: ic

    write(*,'(" [update_temperature] updating temperature")')

    if(count(specific_energy_abs==0.and.density>0.) > 0) then
       write(*,'(" [update_temperature] ",I0," cells have no energy")') count(specific_energy_abs==0.and.density>0.)
    end if

    do ic=1,geo%n_cells
       do id=1,n_dust
          if(is_lte_dust(id)) then
             temperature(ic,id) = specific_energy_abs2temperature(id, specific_energy_abs(ic,id))
          end if
       end do
    end do

    where(temperature < minimum_temperature)
       temperature = minimum_temperature
    end where

    if(allocated(temperature_mean)) call update_mean_temperature()

    call update_energy_abs() ! Update energy absorbed in each cell to match temperature

  end subroutine update_temperature

  subroutine update_energy_abs()
    implicit none
    integer :: ic
    if(main_process()) write(*,'(" [grid_physics] updating energy_abs")')
    do ic=1,geo%n_cells
       do id=1,n_dust
          if(is_lte_dust(id)) then
             specific_energy_abs(ic,id) = temperature2specific_energy_abs(id,temperature(ic,id))
          end if
       end do
    end do
    call update_energy_abs_tot()
  end subroutine update_energy_abs

  subroutine update_energy_abs_tot()
    implicit none
    if(main_process()) write(*,'(" [grid_physics] updating energy_abs_tot")')
    do id=1,n_dust
       energy_abs_tot(id) = sum(specific_energy_abs(:,id)*density(:,id)*geo%volume)
    end do
  end subroutine update_energy_abs_tot

  subroutine precompute_jnu_var()

    implicit none

    integer :: ic,id

    if(main_process()) write(*,'(" [grid_physics] pre-computing jnu_var")')

    do ic=1,geo%n_cells
       do id=1,n_dust   
          call dust_jnu_var_pos_frac(d(id),temperature(ic,id),specific_energy_abs(ic,id),jnu_var_id(ic,id),jnu_var_frac(ic,id))
       end do
    end do

  end subroutine precompute_jnu_var

  real(dp) elemental function difference_ratio(a, b)
    implicit none
    real(dp), intent(in) :: a, b
    difference_ratio = max(a/b, b/a)
  end function difference_ratio

  logical function specific_energy_abs_converged() result(converged)

    implicit none

    real(dp) :: value
    real(dp), save :: value_prev = huge(1._dp)
    real(dp), allocatable, save :: specific_energy_abs_prev(:,:)

    write(*,'(" [specific_energy_abs_converged] checking convergence")')

    if(.not.allocated(specific_energy_abs_prev)) then
       allocate(specific_energy_abs_prev(geo%n_cells, n_dust))
       specific_energy_abs_prev = specific_energy_abs
       converged = .false.
       return
    end if

    value = quantile(reshape(difference_ratio(specific_energy_abs_prev, specific_energy_abs), (/geo%n_cells*n_dust/)), &
         &           convergence_percentile, &
         &           mask=reshape(specific_energy_abs_prev > 0 .and. specific_energy_abs > 0. .and. specific_energy_abs_prev .ne. specific_energy_abs, (/geo%n_cells*n_dust/)))

    write(*,*)
    write(*,'("     -> Percentile: ",F7.2)') convergence_percentile
    write(*,'("     -> Value @ Percentile: ",F)') value
    if(value_prev < huge(1._dp)) then
       write(*,'("     -> Difference from previous iteration: ", F)') difference_ratio(value_prev, value)
    end if
    write(*,*)

    converged = value < convergence_absolute .and. &
         &      abs(difference_ratio(value_prev, value)) < convergence_relative

    specific_energy_abs_prev = specific_energy_abs

    value_prev = value

  end function specific_energy_abs_converged

  type(photon) function emit_from_grid(inu) result(p)

    implicit none

    real(dp) :: xi
    real(dp) :: mass

    integer,intent(in),optional :: inu
    real(dp) :: prob

    if(present(inu)) then

       p%nu = frequencies(inu)
       p%inu = inu

       call prepare_photon(p)
       call update_optconsts(p)

    end if

    call random_number(xi)
    p%dust_id = ceiling(xi*real(n_dust,dp))

    ! Pick random cell
    p%icell = random_cell()            
    p%in_cell = .true.

    ! Find random position inside cell
    call random_position_cell(p%icell, p%r) ! can probably make this a function

    ! Sample an isotropic direction
    call random_sphere_angle3d(p%a)
    call angle3d_to_vector3d(p%a, p%v)

    ! Set stokes parameters to unpolarized light
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    ! Find the relative energy of the photon (relative to grid average energy)
    if(energy_abs_tot(p%dust_id) > 0._dp) then
       mass = density(p%icell%ic,p%dust_id) * geo%volume(p%icell%ic)
       p%energy = specific_energy_abs(p%icell%ic,p%dust_id) &
            & * mass * dble(geo%n_cells) / energy_abs_tot(p%dust_id)
    else
       p%energy = 0._dp
    end if

    ! Set how frequency should be sampled for this photon
    p%emiss_type = 3
    p%emiss_var_id = jnu_var_id(p%icell%ic,p%dust_id)
    p%emiss_var_frac = jnu_var_frac(p%icell%ic,p%dust_id)

    if(present(inu)) then
       call dust_sample_emit_probability(d(p%dust_id),p%emiss_var_id,p%emiss_var_frac,p%nu,prob)
       p%energy = p%energy * prob
    end if

    p%scattered=.false.
    p%reprocessed=.true.
    p%last_isotropic = .true.
    ! p%dust_id already set
    p%last = 'de'

  end function emit_from_grid

end module grid_physics
