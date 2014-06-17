module grid_monochromatic

  ! The purpose of this module is to take care of the emission of photons for
  ! the monochromatic radiative transfer code. The
  ! setup_monochromatic_grid_pdfs subroutine takes care of pre-computing the
  ! probability distribution function for emission at the frequency index inu
  ! (using the frequencies from the radiative transfer settings), and by
  ! default the emission is weighted by the total luminosity from each cell at
  ! that frequency. In future it will be easy to add an option to weigh things
  ! differently, for example by total energy in the cell, or giving equal
  ! weight to each cell. The emit_from_monochromatic_grid_pdf function is then
  ! used to emit a photon from the pre-computed probability distribution
  ! function. The allocate_ and deallocate_ subroutines should be called at
  ! the start and end of the monochormatic radiative transfer respectively.

  use core_lib
  use type_photon, only : photon
  use grid_geometry, only : geo, grid_sample_pdf_map, random_position_cell
  use type_dust, only : dust_sample_emit_probability
  use dust_main, only : d, n_dust, prepare_photon, update_optconsts
  use grid_physics
  use settings, only : frequencies

  implicit none
  save

  private
  public :: allocate_monochromatic_grid_pdfs
  public :: deallocate_monochromatic_grid_pdfs
  public :: setup_monochromatic_grid_pdfs
  public :: emit_from_monochromatic_grid_pdf

  ! Variables for the monochromatic grid emission
  integer :: inu_current
  real(dp), allocatable :: mean_prob(:)
  type(pdf_discrete_dp), allocatable :: emiss_pdf(:)

contains

  subroutine allocate_monochromatic_grid_pdfs()
    implicit none
    allocate(emiss_pdf(n_dust), mean_prob(n_dust))
  end subroutine allocate_monochromatic_grid_pdfs

  subroutine deallocate_monochromatic_grid_pdfs()
    implicit none
    deallocate(emiss_pdf, mean_prob)
  end subroutine deallocate_monochromatic_grid_pdfs

  subroutine setup_monochromatic_grid_pdfs(inu, empty)

    ! Sets up the probability distribution functions to emit from for a specific
    ! frequency for all dust types. This subroutine takes nu, the frequency to
    ! compute the PDF at, and modifies local variables in the module

    implicit none

    integer, intent(in) :: inu

    real(dp) :: nu
    integer :: dust_id
    integer :: emiss_var_id
    real(dp) :: emiss_var_frac
    real(dp), allocatable :: energy(:), prob(:)
    integer :: icell
    logical,intent(out) :: empty

    nu = frequencies(inu)

    ! Allocate temporary arrays
    allocate(energy(geo%n_cells), prob(geo%n_cells))

    ! Loop over dust types
    do dust_id=1, n_dust

       ! Find the total energy emitted inside each cell
       if (energy_abs_tot(dust_id) > 0._dp) then
          energy = specific_energy(:, dust_id) &
               &   * density(:, dust_id) &
               &   * geo%volume(:) &
               &   * dble(geo%n_cells) &
               &   / energy_abs_tot(dust_id)
       else
          energy = 0._dp
       end if

       ! Ensure that energy is zero in masked cells. Note that we can just
       ! work with all the cells here because the masked cells will get
       ! dropped out of the PDF anyway.
       if(geo%masked) then
          where(.not.geo%mask)
             energy = 0._dp
          end where
       end if

       ! Find in each cell the probability of emission at the desired frequency
       ! from the normalized emissivity PDF.
       do icell=1, geo%n_cells
          emiss_var_id = jnu_var_id(icell, dust_id)
          emiss_var_frac = jnu_var_frac(icell, dust_id)
          call dust_sample_emit_probability(d(dust_id), &
               &                            emiss_var_id, emiss_var_frac, &
               &                            nu, prob(icell))
       end do

       ! Set the PDF to the probability of emission at the frequency nu
       mean_prob(dust_id) = mean(prob * energy)
       if(mean_prob(dust_id) > 0._dp) call set_pdf(emiss_pdf(dust_id), prob * energy)

    end do

    ! Deallocate temporary arrays
    deallocate(energy, prob)

    inu_current = inu

    empty = sum(mean_prob) == 0._dp

  end subroutine setup_monochromatic_grid_pdfs

  type(photon) function emit_from_monochromatic_grid_pdf(inu) result(p)

    ! Given a frequency index (inu), sample a position in the grid using the
    ! locally pre-computed PDFs and set up a photon. This function is meant to
    ! be used in conjunction with setup_monochromatic_grid_pdfs, which, given
    ! a frequency, pre-computes the probabilty distribution functon for
    ! emission from the grid.

    implicit none

    integer, intent(in) :: inu

    integer :: dust_id
    real(dp) :: xi

    if(inu /= inu_current) call error('emit_from_monochromatic_grid_pdf', 'incorrect inu')

    p%nu = frequencies(inu)
    p%inu = inu

    call prepare_photon(p)
    call update_optconsts(p)

    call random(xi)
    dust_id = ceiling(xi*real(n_dust, dp))

    if(mean_prob(dust_id) == 0._dp) then
        p%energy = 0._dp
        return
    end if

    call grid_sample_pdf_map(emiss_pdf(dust_id), p%icell)
    p%in_cell = .true.

    ! Find random position inside cell
    call random_position_cell(p%icell, p%r) ! can probably make this a function

    ! Sample an isotropic direction
    call random_sphere_angle3d(p%a)
    call angle3d_to_vector3d(p%a, p%v)

    ! Set stokes parameters to unpolarized light
    p%s = stokes_dp(1._dp, 0._dp, 0._dp, 0._dp)

    ! Set the energy to 1., and it will be scaled in the main routine
    p%energy = mean_prob(dust_id)

    p%scattered = .false.
    p%reprocessed = .true.
    p%last_isotropic = .true.
    p%dust_id = dust_id
    p%last = 'de'

  end function emit_from_monochromatic_grid_pdf

end module grid_monochromatic
