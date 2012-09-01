module type_surface_properties

  use core_lib
  use mpi_hdf5_io
  use type_var2d_pdf2d

  implicit none
  save

  private

  public :: setup_surface_properties

  public :: surface_properties
  type surface_properties

     ! Optical properties

     integer :: n_nu
     real(dp),allocatable :: nu(:)
     real(dp),allocatable :: albedo_nu(:)

     ! Scattering properties
     !
     ! There are four important variables in surface scattering:
     !
     ! * nu: the frequency at which the scattering happens
     ! * i: the incident theta angle (relative to the surface normal)
     ! * e: emerging theta angle (relative to the surface normal)
     ! * psi: the difference in azimuthal angles between the incoming and emergent planes
     !
     ! We use the naming conventions from Hapke (2012). The incoming
     ! photons define nu and mu0, and we then have to sample from a 2-d PDF
     ! that gives the probability of scattering as a function of mu and psi.
     ! We can use the var2d_pdf2d type, which defines a 2-d PDF that
     ! depends on 2 input variables. Bilinear interpolation is used to
     ! interpolate in (nu, mu0), and bilinear interpolation is also used in
     ! the random sampling to interpolate between mu and psi values.

     integer :: n_mu0, n_mu, n_psi
     real(dp),allocatable :: mu0(:), mu(:), psi(:)
     type(var2d_pdf2d_dp) :: radiance

  end type surface_properties

contains

  subroutine setup_surface_properties(group, sp)

    ! Set up surface properties from an HDF5 group
    !
    ! Parameters
    ! ----------
    ! group : HDF5 group/file handle
    !     The group to read the surface properties from
    !
    ! Returns
    ! -------
    ! sp : surface_properties object
    !     The surface properties read from the group

    implicit none

    integer(hid_t), intent(in) :: group
    type(surface_properties), intent(out) :: sp

    real(dp), allocatable :: radiance_array(:, :, :, :)

    character(len=100) :: path

    ! Optical Properties

    path = 'optical_properties'
    call mp_table_read_column_auto(group, path, 'nu', sp%nu)
    call mp_table_read_column_auto(group, path, 'albedo', sp%albedo_nu)
    sp%n_nu = size(sp%nu)

    if(any(is_nan(sp%nu))) call error("setup_surface_properties", "NaN values in nu")

    ! Incident angles

    path = 'incident_angles'
    call mp_table_read_column_auto(group, path, 'mu0', sp%mu0)
    sp%n_mu0 = size(sp%mu0)

    if(any(is_nan(sp%mu0))) call error("setup_surface_properties", "NaN values in i")

    ! Emergent angles

    path = 'emergent_e_angles'
    call mp_table_read_column_auto(group, path, 'mu', sp%mu)
    sp%n_mu = size(sp%mu)

    if(any(is_nan(sp%mu))) call error("setup_surface_properties", "NaN values in e")

    ! Phase angles

    path = 'emergent_psi_angles'
    call mp_table_read_column_auto(group, path, 'psi', sp%psi)
    sp%n_psi = size(sp%psi)

    if(any(is_nan(sp%psi))) call error("setup_surface_properties", "NaN values in psi")

    ! Radiance PDF

    path = 'radiance_pdf'
    call mp_read_array_auto(group, path, radiance_array)

    if(size(radiance_array, 1) /= sp%n_psi) call error("setup_surface_properties", "radiance_array has incorrect dimension 1")
    if(size(radiance_array, 2) /= sp%n_mu) call error("setup_surface_properties", "radiance_array has incorrect dimension 2")
    if(size(radiance_array, 3) /= sp%n_mu0) call error("setup_surface_properties", "radiance_array has incorrect dimension 3")
    if(size(radiance_array, 4) /= sp%n_nu) call error("setup_surface_properties", "radiance_array has incorrect dimension 4")

    if(any(is_nan(radiance_array))) call error("setup_surface_properties", "NaN values in radiance_array")

    sp%radiance = set_var2d_pdf2d(sp%psi, sp%mu, sp%mu0, sp%nu, radiance_array)

  end subroutine setup_surface_properties

end module type_surface_properties
