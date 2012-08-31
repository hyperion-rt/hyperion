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
     ! * g: the scattering phase angle (in the plane of the scattering rays)
     !
     ! We use the naming conventions from Hapke et al. (2009). The incoming 
     ! photons define nu and i, and we then have to sample from a 2-d PDF 
     ! that gives the probability of scattering as a function of e and g. 
     ! We can use the var2d_pdf2d type, which defines a 2-d PDF that 
     ! depends on 2 input variables. Bilinear interpolation is used to 
     ! interpolate in (nu, i), and bilinear interpolation is also used in 
     ! the random sampling to interpolate between e and g values.

     integer :: n_i, n_e, n_g
     real(dp),allocatable :: i(:), e(:), g(:)
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
    call mp_table_read_column_auto(group, path, 'i', sp%i)
    sp%n_i = size(sp%i)

    if(any(is_nan(sp%i))) call error("setup_surface_properties", "NaN values in i")

    ! Emergent angles

    path = 'emergent_angles'
    call mp_table_read_column_auto(group, path, 'e', sp%e)
    sp%n_e = size(sp%e)

    if(any(is_nan(sp%e))) call error("setup_surface_properties", "NaN values in e")

    ! Phase angles

    path = 'phase_angles'
    call mp_table_read_column_auto(group, path, 'g', sp%g)
    sp%n_g = size(sp%g)

    if(any(is_nan(sp%g))) call error("setup_surface_properties", "NaN values in g")

    ! Radiance PDF

    path = 'radiance_pdf'
    call mp_read_array_auto(group, path, radiance_array)

    if(size(radiance_array, 1) /= sp%n_g) call error("setup_surface_properties", "radiance_array has incorrect dimension 1")
    if(size(radiance_array, 2) /= sp%n_e) call error("setup_surface_properties", "radiance_array has incorrect dimension 2")
    if(size(radiance_array, 3) /= sp%n_i) call error("setup_surface_properties", "radiance_array has incorrect dimension 3")
    if(size(radiance_array, 4) /= sp%n_nu) call error("setup_surface_properties", "radiance_array has incorrect dimension 4")

    if(any(is_nan(radiance_array))) call error("setup_surface_properties", "NaN values in radiance_array")

    sp%radiance = set_var2d_pdf2d(sp%g, sp%e, sp%i, sp%nu, radiance_array)

  end subroutine setup_surface_properties

end module type_surface_properties
