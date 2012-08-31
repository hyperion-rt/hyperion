module type_surface

  use core_lib
  use mpi_hdf5_io

  use type_surface_properties, only : surface_properties, setup_surface_properties

  implicit none
  save

  private

  public :: surface_read
  public :: surface_distance
  public :: surface_intersect
  public :: surface_normal

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


end module type_surface
