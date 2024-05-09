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
  use mpi_hdf5_io

  implicit none
  save

  private
  public :: setup_grid_physics
  public :: sublimate_dust
  public :: update_alpha_inv_planck
  public :: check_energy_abs
  public :: update_energy_abs
  public :: update_energy_abs_tot
  public :: select_dust_chi_rho
  public :: select_dust_specific_energy_rho
  public :: emit_from_grid
  public :: precompute_jnu_var
  public :: tau_inv_planck_to_closest_wall
  public :: specific_energy_converged

  ! Density (immutable)
  real(dp),allocatable,target, public :: density(:,:)
  real(dp),allocatable, public :: density_original(:,:)

  ! Variable quantities (made public for MPI)
  integer(idp),allocatable, public :: n_photons(:)
  integer(idp),allocatable, public :: last_photon_id(:)
  real(dp),allocatable, public :: specific_energy(:,:)
  real(dp),allocatable, public :: specific_energy_nu(:,:,:) 
  real(dp),allocatable, public :: specific_energy_sum(:,:)
  real(dp),allocatable, public :: specific_energy_sum_nu(:,:,:)

  real(dp),allocatable, public :: specific_energy_additional(:,:)
  real(dp),allocatable, public :: specific_energy_additional_nu(:,:,:)
  real(dp),allocatable, public :: energy_abs_tot(:)
  real(dp),allocatable, public :: minimum_specific_energy(:)

  real(dp), allocatable,target, public :: alpha_inv_planck(:)

  integer, allocatable, public :: jnu_var_id(:,:)
  real(dp), allocatable, public :: jnu_var_frac(:,:)

  logical, allocatable, public :: empty_neighbor(:)

  ! Temporary variables (used for convenience in external modules)
  real(dp), allocatable, public :: tmp_column_density(:)

  !Counting indices
  integer :: id
  integer :: idx 

  logical :: debug = .false.

  type(pdf_discrete_dp) :: absorption



contains

  real(dp) function tau_inv_planck_to_closest_wall(p) result(tau)
    implicit none
    type(photon),intent(in) :: p
    tau = alpha_inv_planck(p%icell%ic) * distance_to_closest_wall(p)
  end function tau_inv_planck_to_closest_wall

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

  integer function select_dust_specific_energy_rho(icell) result(id_select)
    implicit none
    type(grid_cell),intent(in) :: icell
    do id=1,n_dust
       absorption%pdf(id) = specific_energy(icell%ic, id) * density(icell%ic, id)
    end do
    call find_cdf(absorption)
    id_select = sample_pdf(absorption)
  end function select_dust_specific_energy_rho

  subroutine setup_grid_physics(group, use_mrw, use_pda, compute_isrf)

    implicit none

    integer(hid_t),intent(in) :: group
    logical,intent(in) :: use_mrw, use_pda, compute_isrf
    integer :: n_isrf_wavelengths

    n_isrf_wavelengths = d(1)%n_nu
    ! Density
    allocate(density(geo%n_cells, n_dust))
    allocate(specific_energy(geo%n_cells, n_dust))
    allocate(specific_energy_nu(geo%n_cells,n_dust,n_isrf_wavelengths))

    if(n_dust > 0) then

       if(main_process()) write(*,'(" [grid_physics] reading density grid")')

       ! Read in density
       call read_grid_4d(group, 'density', density, geo)

       ! Check number of dust types for density
       if(size(density, 2).ne.n_dust) call error("setup_grid","density array has wrong number of dust types")

       ! Reset density to zero in masked cells
       if(geo%masked) then
          if(main_process()) write(*, '(" [grid_physics] applying mask to density grid")')
          do id=1,n_dust
             where(.not.geo%mask)
                density(:, id) = 0.
             end where
          end do
       end if

       ! If density difference is requested, save original density
       if(output_density_diff.ne.'none') then
          allocate(density_original(geo%n_cells, n_dust))
          density_original = density
       end if

       if(main_process()) write(*,'(" [grid_physics] reading minimum_specific_energy")')

       ! Read in minimum specific energy
       if(mp_exists_keyword(group, '.', 'minimum_specific_energy')) then
          call mp_read_keyword_vector_auto(group, '.', 'minimum_specific_energy', minimum_specific_energy)
       else
          allocate(minimum_specific_energy(n_dust))
          minimum_specific_energy = 0.
       end if

       if(grid_exists(group, 'specific_energy')) then

          if(main_process()) write(*,'(" [grid_physics] reading specific_energy grid")')

          ! Read in specific_energy
          call read_grid_4d(group, 'specific_energy', specific_energy, geo)

          ! Check number of dust types for specific_energy
          if(size(specific_energy, 2).ne.n_dust) call error("setup_grid","specific_energy array has wrong number of dust types")



          ! Reset specific energy to zero in masked cells
          if(geo%masked) then
             if(main_process()) write(*, '(" [grid_physics] applying mask to specific_energy grid")')
             do id=1,n_dust
                where(.not.geo%mask)
                   specific_energy(:, id) = 0.
                end where
                
                if (compute_isrf) then
                   do idx=1,n_isrf_wavelengths
                      where(.not.geo%mask)
                         specific_energy_nu(:, id, idx) = 0.
                      endwhere
                   end do
                end if

             end do
          end if

          if(trim(specific_energy_type) == 'additional') then
             allocate(specific_energy_additional(geo%n_cells, n_dust))
             if (compute_isrf) then 
                allocate(specific_energy_additional_nu(geo%n_cells,n_dust,n_isrf_wavelengths))
             end if
             ! We store a copy of the initial specific energy in a separate
             ! array, and we set the specific energy to the minimum specific
             ! energy. After the first iteration, specific_energy will get
             ! re-calculated and we will then add specific_energy_additional
             specific_energy_additional = specific_energy
             specific_energy_additional_nu = specific_energy_nu
             do id=1,n_dust
                specific_energy(:,id) = minimum_specific_energy(id)
                
                if (compute_isrf) then
                   do idx=1,n_isrf_wavelengths
                      specific_energy_nu(:,id,idx) = minimum_specific_energy(id)
                   end do
                end if

            end do
          end if


       else

          if(trim(specific_energy_type) == 'additional') then
             call error("setup_grid", "cannot specify specific_energy_type since specific_energy was not given")
          end if

          ! Set all specific_energy to minimum requested
          do id=1,n_dust
             specific_energy(:,id) = minimum_specific_energy(id)

             if (compute_isrf) then
                do idx = 1,n_isrf_wavelengths
                   specific_energy_nu(:,id,idx) = minimum_specific_energy(id)
                end do
             end if

          end do

       end if

    end if

    ! Column density for peeling-off
    allocate(tmp_column_density(n_dust))

    ! Specific energy summation
    allocate(specific_energy_sum(geo%n_cells, n_dust))
    specific_energy_sum = 0._dp

    ! Set up basics for ISRF calculation
    !n_isrf_wavelengths = d(1)%n_nu
    allocate(specific_energy_sum_nu(geo%n_cells, n_dust, n_isrf_wavelengths))
    specific_energy_sum_nu = 0._dp
    
    ! Total energy absorbed
    allocate(energy_abs_tot(n_dust))
    energy_abs_tot = 0._dp

    ! Update energy absorbed in each cell to check bounds and update total
    call check_energy_abs()

    ! Emissivity index and interpolation fraction
    allocate(jnu_var_id(geo%n_cells, n_dust))
    allocate(jnu_var_frac(geo%n_cells, n_dust))

    ! Modified Random Walk
    if(use_mrw) then

       ! Rosseland extinction coefficient
       allocate(alpha_inv_planck(geo%n_cells))
       alpha_inv_planck = 0._dp

    end if

    ! Partial Diffusion Approximation
    if(use_pda .or. trim(output_n_photons) .ne. 'none') then

       ! Number of photons in each cell
       allocate(n_photons(geo%n_cells))
       n_photons = 0

       ! ID of last photon in each cell
       allocate(last_photon_id(geo%n_cells))
       last_photon_id = 0

    end if

    ! Create PDF for absorption in each cell
    call allocate_pdf(absorption,n_dust)

  end subroutine setup_grid_physics

  subroutine update_alpha_inv_planck()

    ! Optimization: could pre-compute alpha_inv_planck just for masked (valid) cells

    implicit none

    integer :: ic

    if(main_process()) write(*,'(" [grid_physics] pre-computing Rosseland absorption coefficient")')
    alpha_inv_planck = 0._dp

    do ic=1,geo%n_cells
       do id=1,n_dust
          if(density(ic, id) > 0._dp) then
             alpha_inv_planck(ic) = alpha_inv_planck(ic) &
                  & + density(ic,id) &
                  & * chi_inv_planck(id, specific_energy(ic,id))
          end if
       end do
    end do

  end subroutine update_alpha_inv_planck

  subroutine sublimate_dust()

    implicit none
    integer :: ic, id
    integer :: reset
    integer :: n_isrf_wavelengths

    n_isrf_wavelengths = d(1)%n_nu

    reset = 0

    do id=1,n_dust

       select case(d(id)%sublimation_mode)
       case(1)

          do ic=1,geo%n_cells
             if(specific_energy(ic, id) > d(id)%sublimation_specific_energy) then
                density(ic, id) = 0.
                specific_energy(ic, id) = minimum_specific_energy(id)
                reset = reset + 1

                if (compute_isrf) then
                   do idx=1,n_isrf_wavelengths
                      specific_energy_nu(ic,id,idx) = minimum_specific_energy(id)
                   end do
                end if

             end if
             
          end do
          if(reset > 0) write(*,'(" [sublimate_dust] dust removed in ",I0," cells")') reset

       case(2)

          do ic=1,geo%n_cells
             if (specific_energy(ic,id) > d(id)%sublimation_specific_energy) then
                density(ic,id) = density(ic,id) &
                     & * d(id)%sublimation_specific_energy / specific_energy(ic, id) &
                     & * (chi_rosseland(id, specific_energy(ic,id)) &
                     & / chi_rosseland(id, d(id)%sublimation_specific_energy))**2
                specific_energy(ic,id) = d(id)%sublimation_specific_energy
                reset = reset + 1


                if (compute_isrf) then
                   do idx=1,n_isrf_wavelengths
                      specific_energy_nu(ic,id,idx) = minimum_specific_energy(id)
                   end do
                end if
                   
             end if
          end do

          if(reset > 0) write(*,'(" [sublimate_dust] density reset due to sublimation in ",I0," cells")') reset

       case(3)

          do ic=1,geo%n_cells
             if(specific_energy(ic, id) > d(id)%sublimation_specific_energy) then
                specific_energy(ic, id) = d(id)%sublimation_specific_energy
                reset = reset + 1
             end if

             
             if (compute_isrf) then
                do idx=1,n_isrf_wavelengths
                   specific_energy_nu(ic,id,idx) = minimum_specific_energy(id)
                end do
             end if

          end do

          if(reset > 0) write(*,'(" [sublimate_dust] capping dust specific_energy in ",I0," cells")') reset

       end select

    end do

    call update_energy_abs_tot()

    call check_energy_abs()

  end subroutine sublimate_dust

  subroutine update_energy_abs(scale)

    implicit none

    real(dp), intent(in) :: scale

    integer :: id,idx
    integer :: n_isrf_wavelengths
    
    n_isrf_wavelengths = d(1)%n_nu

    if(main_process()) write(*,'(" [grid_physics] updating energy_abs")')

    do id=1,n_dust
       specific_energy(:,id) = specific_energy_sum(:,id) * scale / geo%volume
       
       if (compute_isrf) then
          do idx=1,n_isrf_wavelengths
             specific_energy_nu(:,id,idx) = specific_energy_sum_nu(:,id,idx) * scale/geo%volume
          end do
       end if
          
          
       where(geo%volume == 0._dp)
          specific_energy(:,id) = 0._dp
       end where
    end do

    if(count(specific_energy==0.and.density>0.) > 0) then
       write(*,'(" [update_energy_abs] ",I0," cells have no energy")') count(specific_energy==0.and.density>0.)
    end if


    ! Add in additional source of heating
    if(trim(specific_energy_type) == 'additional') then
       if(main_process()) write(*,'(" [grid_physics] adding additional heating source")')
       specific_energy = specific_energy + specific_energy_additional

       
       if (compute_isrf) then 
          specific_energy_nu = specific_energy_nu + specific_energy_additional_nu
       end if

    end if

    call update_energy_abs_tot()

    call check_energy_abs()

  end subroutine update_energy_abs

  subroutine check_energy_abs()

    implicit none

    integer :: id

    if(main_process()) write(*,'(" [grid_physics] checking energy_abs")')

    do id=1,n_dust

       if(any(specific_energy(:,id) < minimum_specific_energy(id))) then
          if(main_process()) call warn("update_energy_abs","specific_energy below minimum requested in some cells - resetting")
          where(specific_energy(:,id) < minimum_specific_energy(id))
             specific_energy(:,id) = minimum_specific_energy(id)
             
 
          end where

       end if

       if(any(specific_energy(:,id) < d(id)%specific_energy(1))) then
          if(enforce_energy_range) then
             if(main_process()) call warn("update_energy_abs","specific_energy below minimum allowed in some cells - resetting")
             where(specific_energy(:,id) < d(id)%specific_energy(1))
                specific_energy(:,id) = d(id)%specific_energy(1)
                
              end where
          else
             if(main_process()) call warn("update_energy_abs","specific_energy below minimum allowed in some cells - will pick closest emissivities")
          end if
       end if

       if(any(specific_energy(:,id) > d(id)%specific_energy(d(id)%n_e))) then
          if(enforce_energy_range) then
             if(main_process()) call warn("update_energy_abs","specific_energy above maximum allowed in some cells - resetting")
             where(specific_energy(:,id) > d(id)%specific_energy(d(id)%n_e))
                specific_energy(:,id) = d(id)%specific_energy(d(id)%n_e)

             end where
          else
             if(main_process()) call warn("update_energy_abs","specific_energy above maximum allowed in some cells - will pick closest emissivities")
          end if
       end if

    end do

    call update_energy_abs_tot()

  end subroutine check_energy_abs

  subroutine update_energy_abs_tot()
    implicit none
    if(main_process()) write(*,'(" [grid_physics] updating energy_abs_tot")')
    do id=1,n_dust
       energy_abs_tot(id) = sum(specific_energy(:,id)*density(:,id)*geo%volume)
    end do
  end subroutine update_energy_abs_tot

  subroutine precompute_jnu_var()

    ! Optimization: could pre-compute jnu_var_id and jnu_var_frac just for masked (valid) cells

    implicit none

    integer :: ic,id

    if(main_process()) write(*,'(" [grid_physics] pre-computing jnu_var")')

    do ic=1,geo%n_cells
       do id=1,n_dust
          call dust_jnu_var_pos_frac(d(id),specific_energy(ic,id),jnu_var_id(ic,id),jnu_var_frac(ic,id))
       end do
    end do

  end subroutine precompute_jnu_var

  real(dp) elemental function difference_ratio(a, b)
    implicit none
    real(dp), intent(in) :: a, b
    difference_ratio = max(a/b, b/a)
  end function difference_ratio

  logical function specific_energy_converged() result(converged)

    implicit none

    real(dp) :: value
    real(dp), save :: value_prev = huge(1._dp)
    real(dp), allocatable, save :: specific_energy_prev(:,:)

    write(*,'(" [specific_energy_converged] checking convergence")')

    if(.not.allocated(specific_energy_prev)) then
       allocate(specific_energy_prev(geo%n_cells, n_dust))
       specific_energy_prev = specific_energy
       converged = .false.
       return
    end if

    if(all(specific_energy_prev .eq. specific_energy)) then
       value = 0._dp
    else if(all(specific_energy_prev .eq. specific_energy .or. specific_energy_prev == 0 .or. specific_energy == 0.)) then
       write(*,*)
       write(*,'(" [specific_energy_converged] could not check for convergence, as the only cells that changed had zero value before or after")')
       converged = .false.
       return
    else
       value = quantile(reshape(difference_ratio(specific_energy_prev, specific_energy), (/geo%n_cells*n_dust/)), &
            &           convergence_percentile, &
            &           mask=reshape(specific_energy_prev > 0 .and. specific_energy > 0. .and. &
            &                        specific_energy_prev .ne. specific_energy, (/geo%n_cells*n_dust/)))
    end if

    write(*,*)
    write(*,'("     -> Percentile: ",F7.2)') convergence_percentile
    write(*,'("     -> Value @ Percentile: ",E10.3)') value
    if(value_prev < huge(1._dp)) then
       if(value == 0._dp) then
          write(*,'("     -> Exact convergence")')
          converged = .true.
       else
          write(*,'("     -> Difference from previous iteration: ", F10.2)') difference_ratio(value_prev, value)
          converged = value < convergence_absolute .and. &
               &      abs(difference_ratio(value_prev, value)) < convergence_relative
       end if
    else
       converged = .false.
    end if
    write(*,*)

    specific_energy_prev = specific_energy

    value_prev = value

  end function specific_energy_converged

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

    call random(xi)
    p%dust_id = ceiling(xi*real(n_dust,dp))

    ! Pick random cell
    p%icell = random_masked_cell()
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
       p%energy = specific_energy(p%icell%ic,p%dust_id) &
            & * mass * dble(geo%n_masked) / energy_abs_tot(p%dust_id)
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
