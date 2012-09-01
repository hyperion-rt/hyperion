module type_surface

  use core_lib
  use mpi_hdf5_io
  use type_var2d_pdf2d
  use type_vector3d
  use type_angle3d
  use type_surface_properties, only : surface_properties, setup_surface_properties

  implicit none
  save

  private

  public :: surface_read
  public :: surface_distance
  public :: surface_intersect
  public :: surface_normal
  public :: surface_scatter
  public :: surface_scatter_peeloff

  public :: surface
  type surface

     ! Surface type. The different surface types implemented so far are:
     !
     ! 2 : sphere

     integer :: type = 0

     ! Surface position, used for spherical surfaces
     type(vector3d_dp) :: position

     ! Surface radius, used for spherical surfaces
     real(dp) :: radius

     ! Surface properties
     type(surface_properties) :: prop

     ! if positive, the surface is connected to a source
     integer :: source_id = -1

  end type surface

contains

  subroutine surface_read(group, srf)

    ! Read a surface from an HDF5 group
    !
    ! Parameters
    ! ----------
    ! group : HDF5 group
    !     The HDF5 group to read the surface from
    !
    ! Returns
    ! -------
    ! srf : surface object
    !     The surface object read in from the HDF5 group

    implicit none

    integer(hid_t),intent(in) :: group
    type(surface), intent(out) :: srf

    integer(hid_t) :: g_prop

    character(len=15) :: type

    ! Read in type first to figure out how to read the rest
    call mp_read_keyword(group, '.', 'type', type)

    select case(trim(type))
    case('sphere')

       srf%type = 2

       call mp_read_keyword(group, '.', 'x', srf%position%x)
       call mp_read_keyword(group, '.', 'y', srf%position%y)
       call mp_read_keyword(group, '.', 'z', srf%position%z)
       call mp_read_keyword(group, '.', 'r', srf%radius)

    case default

       call error("surface_read", "unknown type in surface list: "//trim(type))

    end select

    g_prop = mp_open_group(group, 'surface_properties')
    call setup_surface_properties(g_prop, srf%prop)
    call mp_close_group(g_prop)

  end subroutine surface_read

  real(dp) function surface_distance(srf, r, v)

    ! Find the distance along a ray to a surface
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface to find the distance to
    ! r : vector3d_dp
    !     The current position of the photon
    ! v : vector3d_dp
    !     The direction vector of the photon
    !
    ! Returns
    ! -------
    ! distance : real(dp)
    !     The distance to the surface

    implicit none

    type(surface),intent(in) :: srf
    type(vector3d_dp),intent(in) :: r, v

    type(vector3d_dp) :: dr
    real(dp) :: pB,pC,t1,t2
    real(dp),parameter :: TOL = 1.e-8

    surface_distance = infinity_dp()

    select case(srf%type)
    case(2)

       dr = r - srf%position

       pB = 2._dp * (dr.dot.v)
       pC = (dr.dot.dr) - srf%radius*srf%radius

       call quadratic_pascal_reduced(pB,pC,t1,t2)

       if(t1 < surface_distance .and. t1 > TOL * srf%radius) surface_distance = t1
       if(t2 < surface_distance .and. t2 > TOL * srf%radius) surface_distance = t2

    end select

  end function surface_distance

  logical function surface_intersect(srf, r1, r2)

    ! Find whether a segment (r1, r2) intersects a surface
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface to check for an intersection with
    ! r1, r2 : vector3d_dp
    !     The positions defining the extremities of the segment
    !
    ! Returns
    ! -------
    ! intersect : logical
    !     The distance to the surface

    implicit none

    type(surface), intent(in) :: srf
    type(vector3d_dp), intent(in) :: r1, r2

    type(vector3d_dp) :: r,v
    real(dp) :: A,B,C,t1,t2
    real(dp),parameter :: TOL = 1.e-8

    select case(srf%type)
    case(2)

       r = r1 - srf%position
       v = r2 - r1

       A = v .dot. v
       B = ( r .dot. v ) *  2._dp
       C = ( r .dot. r ) - srf%radius * srf%radius

       call quadratic(A, B, C, t1, t2)

       if((t1 > TOL .and. t1 < 1.) .or. (t2 > TOL .and. t2 < 1.)) then
          surface_intersect = .true.
       else
          surface_intersect = .false.
       end if
    case default
       surface_intersect = .false.
    end select

  end function surface_intersect

  type(vector3d_dp) function surface_normal(srf, r) result(n)

    ! Find the vector normal to a surface
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface to find the normal to
    ! r : vector3d_dp
    !     The position of the photon on the surface
    !
    ! Returns
    ! -------
    ! n : vector3d_dp
    !     The normal to the surface (normalized)

    type(surface), intent(in) :: srf
    type(vector3d_dp), intent(in) :: r
    real(dp) :: norm

    select case(srf%type)
    case(2)

       n = r - srf%position
       n = n / sqrt(n .dot. n)

    end select

  end function surface_normal

  subroutine surface_scatter(srf, nu, r, a, s)

    ! Scatter a photon on a surface
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface the photon is scattering off
    ! nu : real(dp)
    !     The frequency of the incoming photon
    ! r : vector3d_dp
    !     The position of the photon
    ! a : angle3d_dp
    !     The incoming angle
    ! s : stokes object
    !     Incoming Stokes parameters
    !
    ! Returns
    ! -------
    ! a : angle3d_dp
    !     The emergent angle
    ! s : stokes object
    !     Emergent Stokes parameters

    implicit none

    type(surface),intent(in) :: srf
    real(dp),intent(in) :: nu
    type(vector3d_dp),intent(in) :: r
    type(angle3d_dp),intent(inout) :: a
    type(stokes_dp),intent(inout) :: s

    type(angle3d_dp) :: a_local, a_final, a_normal
    type(vector3d_dp) :: v, n
    real(dp) :: mu0, mu, psi

    ! Find local normal vector (normalized)
    n = surface_normal(srf, r)

    ! Convert incoming angle to vector
    call angle3d_to_vector3d(a, v)

    ! Find incident angle
    mu0 = - (v .dot. n)

    ! Check incident angle interval (can't be greater than pi/2)
    if (mu0 < 0._dp .or. mu0 > 1._dp) call error("surface_scatter", "mu0 should be in the range [0:1]")

    ! Sample random outgoing angles
    call sample_var2d_pdf2d(mu0, nu, srf%prop%radiance, psi, mu)

    ! Check emergent angle intervals
    if (psi < 0._dp .or. psi > 2._dp * pi) call error("surface_scatter", "psi should be in the range [0:2pi]")
    if (mu < 0._dp .or. mu > 1._dp) call error("surface_scatter", "mu should be in the range [0:1]")

    ! Construct angle object
    a_local%cosp = cos(psi)
    a_local%sinp = sin(psi)
    a_local%cost = mu
    a_local%sint = sqrt(1._dp - mu * mu)

    call vector3d_to_angle3d(n, a_normal)

    ! Add to original incoming angle
    call rotate_angle3d(a_local, a_normal, a_final)

    ! Set outgoing angle to a_final
    a = a_final

    ! Leave Stokes parameters untouched (for now!)

  end subroutine surface_scatter

  subroutine surface_scatter_peeloff(srf,nu,r,a,s,a_req)

    ! Peeloff a photon from a surface
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface the photon is scattering off
    ! nu : real(dp)
    !     The frequency of the incoming photon
    ! r : vector3d_dp
    !     The position of the photon
    ! a : angle3d_dp
    !     The incoming angle
    ! s : stokes object
    !     Incoming Stokes parameters
    ! a_req : angle3d_dp
    !     The requested peeloff angle
    !
    ! Returns
    ! -------
    ! a : angle3d_dp
    !     The emergent angle
    ! s : stokes object
    !     Emergent Stokes parameters

    implicit none

    type(surface),intent(in) :: srf
    real(dp),intent(in) :: nu
    type(vector3d_dp),intent(in) :: r
    type(angle3d_dp),intent(inout) :: a
    type(stokes_dp),intent(inout) :: s
    type(angle3d_dp),intent(in) :: a_req

    type(angle3d_dp) :: a_local, a_coord
    type(vector3d_dp) :: v, n
    real(dp) :: prob, mu0, psi, mu

    ! Find local normal vector (normalized)
    n = surface_normal(srf, r)

    ! Convert incoming angle to vector
    call angle3d_to_vector3d(a, v)

    ! Recover local angle from normal vector and requested angle
    call vector3d_to_angle3d(n, a_coord)
    call difference_angle3d(a_coord, a_req, a_local)

    if(a_local%cost < 0. .or. a_local%cost > pi / 2.) then

       ! This corresponds to peeloff angles that go through the surface, so
       !the intensities should be set to zero

       s%i = 0.
       s%q = 0.
       s%u = 0.
       s%v = 0.

    else

       ! Find incident angle
       mu0 = - (v .dot. n)

       ! Find emergent e and psi
       mu = a_local%cost
       psi = atan2(a_local%sinp, a_local%cosp)

       ! atan2 gives a value between -pi and pi, so need to wrap
       if(psi < 0._dp) psi = psi + 2. * pi

       ! interpolate the PDF at the requested (i, nu) for the requested (psi,
       ! e). The factor of 4pi is because the rest of the code assumes that
       ! angular distributions are normalized to 4*pi, but the 2-d PDF type
       ! normalizes to 1.

       s%i = interpolate_var2d_pdf2d(mu0, nu, srf%prop%radiance, psi, mu) * 4._dp * pi
       s%q = 0.
       s%u = 0.
       s%v = 0.

    end if

    a = a_req

  end subroutine surface_scatter_peeloff

end module type_surface
