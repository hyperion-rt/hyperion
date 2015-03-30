module type_source

  use core_lib
  use mpi_hdf5_io
  use type_grid_cell
  use type_photon
  use grid_geometry
  use grid_physics
  use type_dust
  use dust_main, only : d

  implicit none
  save

  private
  public :: source
  public :: source_read
  public :: source_emit
  public :: source_emit_peeloff
  public :: source_intersect
  public :: source_distance
  public :: get_spectrum_interp
  public :: get_spectrum_binned

  public :: spot

  type spot
     type(angle3d_dp) :: a
     real(dp) :: cost
     integer :: freq_type = 0
     type(pdf_dp) :: spectrum
     real(dp) :: temperature
  end type spot

  type source

     ! Source type. The different source types implemented so far are:
     !
     ! 1 : point source
     ! 2 : sphere
     ! 3 : sphere section (spot)
     ! 4 : luminosity grid (used e.g. for accretion disk luminosity)
     ! 5 : external illumination sphere
     ! 6 : external illumination box
     ! 7 : plane parallel source
     ! 8 : point source collection

     integer  :: type = 0

     ! Source luminosity
     real(dp) :: luminosity

     ! Whether the source should be peeled off
     logical :: peeloff

     ! Source position, used for point, spherical, and spot sources
     type(vector3d_dp) :: position

     ! Collection of positions (used for PointSourceCollection)
     type(pdf_discrete_dp) :: collection_pdf
     type(vector3d_dp),allocatable :: position_collection(:)

     ! Spot position and size
     integer :: n_spots = 0
     type(pdf_discrete_dp) :: spot_pdf
     type(spot),allocatable :: spot(:)

     ! Source radius and limb darkening, used for spherical and spot sources
     real(dp) :: radius
     logical  :: limb_darkening

     ! Range for external box source
     real(dp) :: xmin, xmax
     real(dp) :: ymin, ymax
     real(dp) :: zmin, zmax
     type(pdf_discrete_dp) :: face

     ! For plane parallel
     type(angle3d_dp) :: direction

     ! Luminosity grid
     type(pdf_discrete_dp) :: luminosity_map

     ! How to sample frequencies. The different types implemented so far are:
     !
     ! 1 : spectrum
     ! 2 : blackbody
     ! 3 : local emissivity

     integer :: freq_type = 0

     ! Source spectrum, used for spectral sampling
     type(pdf_dp) :: spectrum
     real(dp) :: temperature

     logical :: intersect = .false.

  end type source

contains

  subroutine source_read(group, s)
    implicit none
    integer(hid_t),intent(in) :: group
    type(source), intent(out) :: s
    character(len=50) :: type
    real(dp) :: lon, lat, spot_size, luminosity
    integer :: i
    character(len=255),allocatable :: spot_names(:)
    integer(hid_t) :: g_spot
    real(dp) :: dx, dy, dz
    real(dp) :: theta, phi
    real(dp),allocatable :: luminosity_collection(:), position_collection(:,:)

    call mp_read_keyword(group, '.', 'type', type)

    select case(trim(type))
    case('point')

       s%type = 1

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_keyword(group, '.', 'x', s%position%x)
       call mp_read_keyword(group, '.', 'y', s%position%y)
       call mp_read_keyword(group, '.', 'z', s%position%z)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "Point source cannot have LTE spectrum")

    case('sphere')

       s%type = 2

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_keyword(group, '.', 'x', s%position%x)
       call mp_read_keyword(group, '.', 'y', s%position%y)
       call mp_read_keyword(group, '.', 'z', s%position%z)
       call mp_read_keyword(group, '.', 'r', s%radius)
       call mp_read_keyword(group, '.', 'limb', s%limb_darkening)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "Spherical source cannot have LTE spectrum")

       s%intersect = .true.

       call mp_list_groups(group, '.', spot_names)
       s%n_spots = size(spot_names)

       ! If there are no spots, we are done
       if(s%n_spots==0) return

       ! Otherwise change source type to spotted sphere
       s%type = 3

       ! Allocate PDF for sampling from spots/sphere
       call allocate_pdf(s%spot_pdf,s%n_spots+1)
       s%spot_pdf%pdf(s%n_spots + 1) = s%luminosity

       ! Allocate arrays for spot angle and size
       allocate(s%spot(s%n_spots))

       ! Read in the spot parameters
       do i=1,s%n_spots

          call mp_read_keyword(group, spot_names(i), 'luminosity', luminosity)
          call mp_read_keyword(group, spot_names(i), 'longitude', lon)
          call mp_read_keyword(group, spot_names(i), 'latitude', lat)
          call mp_read_keyword(group, spot_names(i), 'radius', spot_size)

          s%spot_pdf%pdf(i) = luminosity
          s%luminosity = s%luminosity + luminosity
          s%spot(i)%a = angle3d_deg(lon, lat)
          s%spot(i)%cost = cos(spot_size * deg2rad)

          g_spot = mp_open_group(group, spot_names(i))
          call set_spectrum(group, s%spot(i)%freq_type, s%spot(i)%spectrum, s%spot(i)%temperature)
          call mp_close_group(g_spot)

          if(s%freq_type == 3) call error("source_read", "Spot cannot have LTE spectrum")

       end do

       ! Compute CDF for spots/sphere selection PDF
       call find_cdf(s%spot_pdf)

    case('map')

       s%type = 4

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       call grid_load_pdf_map(group, 'Luminosity map', s%luminosity_map)

    case('extern_sph')

       s%type = 5

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_keyword(group, '.', 'x', s%position%x)
       call mp_read_keyword(group, '.', 'y', s%position%y)
       call mp_read_keyword(group, '.', 'z', s%position%z)
       call mp_read_keyword(group, '.', 'r', s%radius)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "External spherical source cannot have LTE spectrum")

    case('extern_box')

       s%type = 6

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_keyword(group, '.', 'xmin', s%xmin)
       call mp_read_keyword(group, '.', 'xmax', s%xmax)
       call mp_read_keyword(group, '.', 'ymin', s%ymin)
       call mp_read_keyword(group, '.', 'ymax', s%ymax)
       call mp_read_keyword(group, '.', 'zmin', s%zmin)
       call mp_read_keyword(group, '.', 'zmax', s%zmax)

       dx = s%xmax - s%xmin
       dy = s%ymax - s%ymin
       dz = s%zmax - s%zmin

       call set_pdf(s%face, (/dy*dz, dy*dz, dz*dx, dz*dx, dx*dy, dx*dy/))

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "External box source cannot have LTE spectrum")

    case('plane_parallel')

       s%type = 7

       call mp_read_keyword(group, '.', 'luminosity', s%luminosity)
       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_keyword(group, '.', 'x', s%position%x)
       call mp_read_keyword(group, '.', 'y', s%position%y)
       call mp_read_keyword(group, '.', 'z', s%position%z)
       call mp_read_keyword(group, '.', 'r', s%radius)
       call mp_read_keyword(group, '.', 'theta', theta)
       call mp_read_keyword(group, '.', 'phi', phi)

       s%direction = angle3d_deg(theta, phi)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "Plane parallel cannot have LTE spectrum")

    case('point_collection')

       s%type = 8

       call mp_read_keyword(group, '.', 'peeloff', s%peeloff)
       call mp_read_array_auto(group, 'position', position_collection)
       call mp_read_array_auto(group, 'luminosity', luminosity_collection)

       allocate(s%position_collection(size(position_collection, 2)))

       s%position_collection%x = position_collection(1,:)
       s%position_collection%y = position_collection(2,:)
       s%position_collection%z = position_collection(3,:)

       s%luminosity = sum(luminosity_collection)
       call set_pdf(s%collection_pdf, luminosity_collection)

       call set_spectrum(group, s%freq_type, s%spectrum, s%temperature)

       if(s%freq_type == 3) call error("source_read", "Point source collection cannot have LTE spectrum")

    case default
       call error("source_read", "unknown type in source list: "//trim(type))
    end select
  end subroutine source_read


  subroutine set_spectrum(group, freq_type, spectrum, temperature)

    implicit none

    integer(hid_t),intent(in) :: group

    type(pdf_dp),intent(out) :: spectrum
    real(dp),intent(out) :: temperature
    integer,intent(out) :: freq_type

    character(len=255) :: spec_type

    real(dp),allocatable :: nu(:), fnu(:)
    integer :: inu

    call mp_read_keyword(group, '.', 'spectrum', spec_type)

    select case(trim(spec_type))
    case('spectrum')
       call mp_table_read_column_auto(group, 'spectrum', 'nu', nu)
       call mp_table_read_column_auto(group, 'spectrum', 'fnu', fnu)
       do inu=1,size(nu)-1
          if(nu(inu + 1) < nu(inu)) then
             call error("set_spectrum", "spectrum frequency should be monotonically increasing")
          end if
       end do
       call set_pdf(spectrum,nu,fnu,log=.true.)
       freq_type = 1
    case('temperature')
       call mp_read_keyword(group, '.', 'temperature', temperature)
       freq_type = 2
    case('lte')
       freq_type = 3
    case default
       call error("set_spectrum", "unknown spectrum specifier: "//trim(spec_type))
    end select

  end subroutine set_spectrum

  real(dp) function source_distance(src,r,v)

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    type(vector3d_dp),intent(in) :: r,v
    ! photon position and direction

    type(vector3d_dp) :: dr

    real(dp) :: pB,pC,t1,t2

    real(dp),parameter :: tol = 1.e-8

    source_distance = infinity_dp()

    select case(src%type)
    case(2,3)

       dr = r - src%position

       pB = 2._dp * (dr.dot.v)
       pC = (dr.dot.dr) - src%radius*src%radius

       call quadratic_pascal_reduced(pB,pC,t1,t2)

       if(t1 < source_distance.and.t1 > tol*src%radius) source_distance = t1
       if(t2 < source_distance.and.t2 > tol*src%radius) source_distance = t2

    end select

  end function source_distance

  logical function source_intersect(src,r1,r2)

    implicit none
    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    type(vector3d_dp),intent(in) :: r1,r2
    ! start and end position

    type(vector3d_dp) :: r,v

    real(dp) :: A,B,C,t1,t2

    select case(src%type)
    case(2)

       r = r1 - src%position
       v = r2 - r1

       A = v .dot. v
       B = ( r .dot. v )* 2._dp
       C = ( r .dot. r ) - src%radius*src%radius

       call quadratic(A,B,C,t1,t2)

       if((t1.gt.1.e-8.and.t1.lt.1.).or.&
            &(t2.gt.1.e-8.and.t2.lt.1.)) then
          source_intersect = .true.
       else
          source_intersect = .false.
       end if
    case default
       source_intersect = .false.
    end select

  end function source_intersect

  subroutine source_emit(src,p,nu)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    real(dp),intent(in),optional :: nu
    ! frequency required - if this is set, then the stokes intensity reflects the probability of emission at this wavelength

    ! --- Output --- !
    type(photon),intent(inout) :: p
    ! the position and direction of the emitted photon

    integer :: ispot

    select case(src%type)
    case(1)
       call emit_from_point(src,p)
    case(2)
       call emit_from_sphere(src,p)
    case(3)
       ispot = sample_pdf(src%spot_pdf)
       if(ispot==src%n_spots+1) then
          call emit_from_sphere(src,p)
       else
          call emit_from_sphere(src,p,spot=ispot)
       end if
       ! Need to change ispot if sphere and falls inside spot
    case(4)
       call emit_from_map(src,p)
    case(5)
       call emit_from_extern_sph(src,p)
    case(6)
       call emit_from_extern_box(src,p)
    case(7)
       call emit_from_plane_parallel(src,p)
    case(8)
       call emit_from_point_collection(src,p)
    end select

    if(present(nu)) then

       ! The frequency is fixed, but the Stokes intensity needs to be modified to reflect the probability of emission

       p%nu = nu

       if(src%type==3.and.ispot .le. src%n_spots) then
          select case(src%spot(ispot)%freq_type)
          case(1)
             p%energy = interpolate_pdf(src%spot(ispot)%spectrum, nu, bounds_error=.false., fill_value=0._dp)
          case(2)
             p%energy = normalized_B_nu(nu, src%spot(ispot)%temperature)
          case(3)
             p%dust_id = select_dust_specific_energy_rho(p%icell)
             p%emiss_var_id = jnu_var_id(p%icell%ic, p%dust_id)
             p%emiss_var_frac = jnu_var_frac(p%icell%ic, p%dust_id)
             call dust_sample_emit_probability(d(p%dust_id),p%emiss_var_id,p%emiss_var_frac,nu,p%energy)
          end select
       end if

       select case(src%freq_type)
       case(1)
          p%energy = interpolate_pdf(src%spectrum, nu, bounds_error=.false., fill_value=0._dp)
       case(2)
          p%energy = normalized_B_nu(nu, src%temperature)
       case(3)
          p%dust_id = select_dust_specific_energy_rho(p%icell)
          p%emiss_var_id = jnu_var_id(p%icell%ic, p%dust_id)
          p%emiss_var_frac = jnu_var_frac(p%icell%ic, p%dust_id)
          call dust_sample_emit_probability(d(p%dust_id),p%emiss_var_id,p%emiss_var_frac,nu,p%energy)
       end select

    else

       p%energy = 1._dp

       if(src%type==3.and.ispot .le. src%n_spots) then
          select case(src%spot(ispot)%freq_type)
          case(1)
             p%nu = sample_pdf_log(src%spot(ispot)%spectrum)
          case(2)
             call random_planck_frequency(p%nu, src%spot(ispot)%temperature)
          case(3)
             p%dust_id = select_dust_specific_energy_rho(p%icell)
             p%emiss_var_id = jnu_var_id(p%icell%ic, p%dust_id)
             p%emiss_var_frac = jnu_var_frac(p%icell%ic, p%dust_id)
             call dust_sample_j_nu(d(p%dust_id),p%emiss_var_id,p%emiss_var_frac,p%nu)
          end select
       end if

       select case(src%freq_type)
       case(1)
          p%nu = sample_pdf(src%spectrum)
       case(2)
          call random_planck_frequency(p%nu, src%temperature)
       case(3)
          p%dust_id = select_dust_specific_energy_rho(p%icell)
          p%emiss_var_id = jnu_var_id(p%icell%ic, p%dust_id)
          p%emiss_var_frac = jnu_var_frac(p%icell%ic, p%dust_id)
          call dust_sample_j_nu(d(p%dust_id),p%emiss_var_id,p%emiss_var_frac,p%nu)
       end select

    end if

  end subroutine source_emit

  subroutine source_emit_peeloff(src,p,a_req)
    implicit none
    type(source),intent(in) :: src ! the source to emit from
    type(photon),intent(inout) :: p ! the photon to peeloff
    type(angle3d_dp),intent(in) :: a_req ! requested angle
    if(src%peeloff) then
       select case(src%type)
       case(2,3)
          call emit_from_sphere_peeloff(src,p,a_req)
       case(5)
          call emit_from_extern_sph_peeloff(src,p,a_req)
       case(6)
          call emit_from_extern_box_peeloff(src,p,a_req)
       case default
          stop "Should not be here, all other source types are isotropic"
       end select
    else
       p%s = stokes_dp(0._dp, 0._dp, 0._dp, 0._dp)
       p%a = a_req
    end if
  end subroutine source_emit_peeloff

  !**********************************************************************!
  ! emit_from_point : emit a photon from a point source s
  !**********************************************************************!

  subroutine emit_from_point(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    ! Set position to that of point source
    p%r = src%position

    ! Sample isotropic angle
    call random_sphere_angle3d(p%a)

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    p%last_isotropic = .true.

  end subroutine emit_from_point

  !**********************************************************************!
  ! emit_from_point_collection : emit a photon from a point source collection
  !**********************************************************************!

  subroutine emit_from_point_collection(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    integer :: i_source

    ! Set position to that of one of the point sources
    i_source = sample_pdf(src%collection_pdf)
    p%r = src%position_collection(i_source)

    ! Sample isotropic angle
    call random_sphere_angle3d(p%a)

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    p%last_isotropic = .true.

  end subroutine emit_from_point_collection

  !**********************************************************************!
  ! emit_from_sphere : emit a photon from a sphere s
  !**********************************************************************!

  subroutine emit_from_sphere(src,p,spot)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    integer,intent(in),optional :: spot

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    ! --- Local variables --- !

    type(angle3d_dp)  :: a_coord,a_local

    ! temporary position

    real(dp) :: phi_local

    real(dp) :: xi

    ! --- First, generate a random position on a unit sphere --- !

    if(present(spot)) then
       do
          call random_sphere_angle3d(a_coord)
          if((a_coord.dot.src%spot(spot)%a) .gt. src%spot(spot)%cost) exit
       end do
    else
       call random_sphere_angle3d(a_coord)
    end if

    ! --- Sample a random direction for photon emission --- !
    !
    ! Because we are emitting from the surface of a star,
    ! we need to sample the (local) phi and theta angles
    ! from:
    !
    ! P(phi)   = 1/(2pi) * dphi
    ! P(theta) = 2 * mu * dmu    where mu = cos(theta)
    !
    ! This means we have:
    !
    ! phi   = xi * 2*pi
    ! cos_theta = sqrt(xi)
    !
    ! where xi is a random number between 0 and 1

    call random_uni(phi_local,zero,twopi)

    a_local%cosp = cos(phi_local)
    a_local%sinp = sin(phi_local)

    if(src%limb_darkening) then
       a_local%cost = ran_mu_limb(1.5_dp,1.0_dp)
    else
       call random(xi)
       a_local%cost = sqrt(xi)
    end if

    a_local%sint = sqrt( 1._dp - a_local%cost * a_local%cost )

    ! --- And rotate to general frame of reference --- !

    call rotate_angle3d(a_local,a_coord,p%a)

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    ! --- Convert position angle on star to real position --- !

    call angle3d_to_vector3d(a_coord,p%r)

    ! --- Set positon to stellar position + vector to surface --- !

    p%r = p%r * src%radius
    p%r = p%r + src%position

    p%last_isotropic = .false.
    p%source_a = a_coord

  end subroutine emit_from_sphere

  subroutine emit_from_sphere_peeloff(src,p,a_req)
    implicit none
    type(source),intent(in) :: src ! the source to emit from
    type(photon),intent(inout) :: p ! the photon to peeloff
    type(angle3d_dp),intent(in) :: a_req ! requested angle
    real(dp) :: mu
    mu = max(a_req .dot. p%source_a, 0._dp)
    ! The probability distribution functions are normalized so that their
    ! total integrals are 4*pi (not 1)
    if(src%limb_darkening) then
       p%s = stokes_dp(2.*(1.5_dp * mu*mu + mu), 0._dp, 0._dp, 0._dp)
    else
       p%s = stokes_dp(4.*mu, 0._dp, 0._dp, 0._dp)
    end if
    p%a = a_req
  end subroutine emit_from_sphere_peeloff

  !**********************************************************************!
  ! emit_from_map : emit a photon from a luminosity map
  !**********************************************************************!

  subroutine emit_from_map(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    ! Sample position in map
    call grid_sample_pdf_map(src%luminosity_map, p%icell)
    p%in_cell = .true.

    ! Find random position inside cell
    call random_position_cell(p%icell, p%r)

    ! Sample isotropic angle
    call random_sphere_angle3d(p%a)

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    p%last_isotropic = .true.

  end subroutine emit_from_map

  !**********************************************************************!
  ! emit_from_extern_sph : emit a photon from a sphere s, inwards
  !**********************************************************************!

  subroutine emit_from_extern_sph(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    ! --- Local variables --- !

    type(angle3d_dp)  :: a_coord,a_local

    ! temporary position

    real(dp) :: phi_local

    real(dp) :: xi

    ! --- First, generate a random position on a unit sphere --- !

    call random_sphere_angle3d(a_coord)

    ! See emit_sphere for explanation

    call random_uni(phi_local,zero,twopi)

    a_local%cosp = cos(phi_local)
    a_local%sinp = sin(phi_local)

    call random(xi)
    a_local%cost = sqrt(xi)
    a_local%sint = sqrt( 1._dp - a_local%cost * a_local%cost )

    ! --- And rotate to general frame of reference --- !

    call rotate_angle3d(a_local,a_coord,p%a)

    ! Point inwards
    p%a = - p%a

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    ! --- Convert position angle on sphere to real position --- !

    call angle3d_to_vector3d(a_coord,p%r)

    ! --- Set positon to stellar position + vector to surface --- !

    p%r = p%r * src%radius
    p%r = p%r + src%position

    p%last_isotropic = .false.
    p%source_a = -a_coord

  end subroutine emit_from_extern_sph

  subroutine emit_from_extern_sph_peeloff(src,p,a_req)
    implicit none
    type(source),intent(in) :: src ! the source to emit from
    type(photon),intent(inout) :: p ! the photon to peeloff
    type(angle3d_dp),intent(in) :: a_req ! requested angle
    real(dp) :: mu
    mu = max(a_req .dot. p%source_a, 0._dp)
    p%s = stokes_dp(4._dp*mu, 0._dp, 0._dp, 0._dp)
    p%a = a_req
  end subroutine emit_from_extern_sph_peeloff

  subroutine emit_from_extern_box(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    ! --- Local variables --- !

    type(angle3d_dp)  :: a_coord,a_local

    ! temporary position

    real(dp) :: phi_local

    real(dp) :: xi

    integer :: face

    ! Sample a side to emit from
    face = sample_pdf(src%face)

    ! Sample direction vector

    call random_uni(phi_local,zero,twopi)
    a_local%cosp = cos(phi_local)
    a_local%sinp = sin(phi_local)

    call random(xi)
    a_local%cost = sqrt(xi)
    a_local%sint = sqrt( 1._dp - a_local%cost * a_local%cost )

    ! Set position and rotation angle depending on the face

    select case(face)
    case(1)
       p%r%x = src%xmin
       call random_uni(p%r%y, src%ymin, src%ymax)
       call random_uni(p%r%z, src%zmin, src%zmax)
       a_coord = angle3d_dp(0._dp, 1._dp, 1._dp, 0._dp)
    case(2)
       p%r%x = src%xmax
       call random_uni(p%r%y, src%ymin, src%ymax)
       call random_uni(p%r%z, src%zmin, src%zmax)
       a_coord = angle3d_dp(0._dp, -1._dp, 1._dp, 0._dp)
    case(3)
       call random_uni(p%r%x, src%xmin, src%xmax)
       p%r%y = src%ymin
       call random_uni(p%r%z, src%zmin, src%zmax)
       a_coord = angle3d_dp(0._dp, 1._dp, 0._dp, 1._dp)
    case(4)
       call random_uni(p%r%x, src%xmin, src%xmax)
       p%r%y = src%ymax
       call random_uni(p%r%z, src%zmin, src%zmax)
       a_coord = angle3d_dp(0._dp, -1._dp, 0._dp, 1._dp)
    case(5)
       call random_uni(p%r%x, src%xmin, src%xmax)
       call random_uni(p%r%y, src%ymin, src%ymax)
       p%r%z = src%zmin
       a_coord = angle3d_dp(1._dp, 0._dp, 1._dp, 0._dp)
    case(6)
       call random_uni(p%r%x, src%xmin, src%xmax)
       call random_uni(p%r%y, src%ymin, src%ymax)
       p%r%z = src%zmax
       a_coord = angle3d_dp(-1._dp, 0._dp, 1._dp, 0._dp)
    end select

    call rotate_angle3d(a_local,a_coord,p%a)

    ! Set Stokes vector
    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    ! --- Set positon to stellar position + vector to surface --- !

    p%last_isotropic = .false.

    p%face_id = face

  end subroutine emit_from_extern_box

  subroutine emit_from_extern_box_peeloff(src,p,a_req)
    implicit none
    type(source),intent(in) :: src ! the source to emit from
    type(photon),intent(inout) :: p ! the photon to peeloff
    type(angle3d_dp),intent(in) :: a_req ! requested angle
    real(dp) :: mu
    type(angle3d_dp) :: a_coord
    select case(p%face_id)
    case(1)
       a_coord = angle3d_dp(0._dp, 1._dp, 1._dp, 0._dp)
    case(2)
       a_coord = angle3d_dp(0._dp, -1._dp, 1._dp, 0._dp)
    case(3)
       a_coord = angle3d_dp(0._dp, 1._dp, 0._dp, 1._dp)
    case(4)
       a_coord = angle3d_dp(0._dp, -1._dp, 0._dp, 1._dp)
    case(5)
       a_coord = angle3d_dp(1._dp, 0._dp, 1._dp, 0._dp)
    case(6)
       a_coord = angle3d_dp(-1._dp, 0._dp, 1._dp, 0._dp)
    end select
    mu = max(a_req .dot. a_coord, 0._dp)
    p%s = stokes_dp(4._dp*mu, 0._dp, 0._dp, 0._dp)
    p%a = a_req
  end subroutine emit_from_extern_box_peeloff

  subroutine emit_from_plane_parallel(src,p)

    implicit none

    ! --- Input --- !

    type(source),intent(in) :: src
    ! the source to emit from

    ! --- Output --- !

    type(photon),intent(inout) :: p
    ! the emitted photon

    real(dp) :: xi, r, phi

    type(angle3d_dp) :: a_local, a_final

    ! Sample radius on emitting disk
    call random(xi)
    r = xi**0.5 * src%radius

    ! Sample phi on emitting disk
    call random_uni(phi, 0._dp, 360._dp)

    ! Find coordinate angle of point on disk
    a_local = angle3d_deg(90._dp, phi)
    call rotate_angle3d(a_local,src%direction,a_final)

    ! Convert to position, and apply radius
    call angle3d_to_vector3d(a_final, p%r)
    p%r = p%r * r
    p%r = p%r + src%position

    ! Set photon direction
    p%a = src%direction

    p%s = stokes_dp(1._dp,0._dp,0._dp,0._dp)

    p%last_isotropic = .false.

  end subroutine emit_from_plane_parallel

  !**********************************************************************!
  ! ran_mu_limb(a,b) : sample random mu with limb darkening a*mu^2+b*mu
  !**********************************************************************!

  real(dp) function ran_mu_limb(a,b)

    implicit none

    ! Returns a random number sampled from the limb
    ! darkened distribution with PDF
    !
    ! P = a*mu**2 + b*mu

    real(dp),intent(in) :: a,b
    ! the coefficients of the probability distribution

    real(dp) :: s,t,xi,norm
    ! temporary variables

    real(dp),parameter :: half  = 1._dp/2._dp
    real(dp),parameter :: third = 1._dp/3._dp

    ! The CDF is
    ! CDF = a/3 * mu**3 + b/2 * mu**2 = s*mu**3+t*mu**2

    s = a * third
    t = b * half

    norm = s + t

    s = s / norm
    t = t / norm

    ! The equation to solve is xi = s*mu**3+t*mu**2

    ! So we want to solve the cubic equation:
    !
    ! s*mu**3 + t*mu**2 - xi = 0
    !
    ! Or
    !
    ! mu**3 + t/s * mu**2 - xi/s = 0

    call random(xi)
    xi = -xi

    ran_mu_limb = cubic_real_root_v2(t/s,xi/s)

  contains

    real(dp) function cubic_real_root_v2(b,d)

      implicit none

      ! solves the quadratic equation
      !
      ! x^3 + b*x**2 + d = 0

      real(dp) :: b,d
      real(dp) :: p,q,p3,q2

      real(dp),parameter :: alpha = 1._dp/3._dp
      real(dp),parameter :: beta  = 1._dp/9._dp
      real(dp),parameter :: gamma = 1._dp/27._dp

      real(dp) :: delta
      real(dp) :: u,v,y

      real(dp) :: phi

      ! Compute p and q
      ! Note - here I define q as the usual q divided by 2

      p = - b*b*alpha* alpha
      q = ( d + 2._dp*b*b*b*gamma ) * 0.5_dp

      p3 = p*p*p
      q2 = q*q

      ! Compute discriminant

      delta = q2 + p3

      if(delta < 0) then

         phi = acos(-q/sqrt(abs(p3)))

         y = + 2 * sqrt(abs(p)) * cos(phi*alpha)

         cubic_real_root_v2 = y - b*alpha

      else

         ! Compute u and v

         delta = sqrt(delta)

         u = cbrt(-q+delta)
         v = cbrt(-q-delta)

         ! Find real root

         cubic_real_root_v2 = u + v - b*alpha

      end if

    end function cubic_real_root_v2

  end function ran_mu_limb

  elemental real(dp) function normalized_B_nu(nu,T)
    implicit none
    real(dp),intent(in) :: nu,T
    real(dp),parameter :: a = two * h_cgs / c_cgs / c_cgs / stef_boltz * pi
    real(dp),parameter :: b = h_cgs / k_cgs
    real(dp) :: T4
    T4 = T*T*T*T
    normalized_B_nu = a * nu * nu * nu / ( exp(b*nu/T) - one) / T4
  end function normalized_B_nu

  function get_spectrum_interp(src, nu) result(spectrum)

    implicit none

    type(source),intent(in) :: src
    real(dp),intent(in) :: nu(:)
    real(dp) :: spectrum(size(nu))

    select case(src%freq_type)
    case(1)
       spectrum = interp1d_loglog(src%spectrum%x, src%spectrum%pdf, &
            &                 nu, bounds_error=.false., fill_value=0._dp)
    case(2)
       spectrum = normalized_B_nu(nu, src%temperature)
    case default
       call error("get_spectrum_interp", "cannot get spectrum")
    end select

  end function get_spectrum_interp

  function get_spectrum_binned(src, n_nu, nu_min, nu_max) result(spectrum)

    implicit none

    type(source),intent(in) :: src
    integer,intent(in) :: n_nu
    real(dp),intent(in) :: nu_min, nu_max
    real(dp) :: spectrum(n_nu)

    real(dp),allocatable :: nu(:), fnu(:)
    integer :: inu, n_nu_bb
    real(dp) :: numin, numax
    real(dp) :: log10_nu_min_bb, log10_nu_max_bb
    real(dp) :: log10_nu_min, log10_nu_max

    select case(src%freq_type)
    case(1)

       allocate(nu(src%spectrum%n))
       allocate(fnu(src%spectrum%n))

       nu = src%spectrum%x
       fnu = src%spectrum%pdf

    case(2)

       log10_nu_min_bb = log10(3.e9_dp)  ! 10 cm (0.05K)
       log10_nu_max_bb = log10(3.e16_dp) ! 10 nm (1/2 million K)

       n_nu_bb = ceiling((log10_nu_max_bb - log10_nu_min_bb) * 100000)
       allocate(nu(n_nu_bb))
       allocate(fnu(n_nu_bb))
       do inu=1,n_nu_bb
          nu(inu) = 10._dp ** (real(inu-1,dp) / real(n_nu_bb-1,dp)*(log10_nu_max_bb - log10_nu_min_bb) + log10_nu_min_bb)
       end do

       fnu = normalized_B_nu(nu, src%temperature)

    end select

    log10_nu_min = log10(nu_min)
    log10_nu_max = log10(nu_max)

    do inu=1, n_nu

       numin = 10._dp ** (log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu-1, dp) / real(n_nu, dp))
       numax = 10._dp ** (log10_nu_min + (log10_nu_max - log10_nu_min) * real(inu, dp) / real(n_nu, dp))

       spectrum(inu) = integral_loglog(nu, fnu, numin, numax)

    end do

    spectrum = spectrum / integral_loglog(nu, fnu)

  end function get_spectrum_binned

end module type_source
