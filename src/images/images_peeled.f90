module peeled_images

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use grid_geometry, only : place_in_cell
  use grid_propagate

  use type_image
  use type_photon
  use type_source
  use type_dust

  use sources
  use dust_main
  use dust_interact

  use grid_mrw

  implicit none
  save

  private

  public :: peeled_images_adjust_scale
  public :: peeled_images_setup
  public :: peeled_images_write
  public :: peeloff_photon
  public :: a_peeloff

  logical,public :: make_peeled_images

  integer,public :: n_groups = 0
  integer,public :: n_peeled = 0
  ! number of directions

  integer,allocatable :: group_id(:)
  integer,allocatable :: view_id(:)

  logical,allocatable :: inside_observer(:)
  logical,allocatable :: ignore_optical_depth(:)
  type(vector3d_dp), allocatable :: r_peeloff(:)
  type(angle3d_dp), allocatable :: viewing_angles(:)
  ! angles in which to make peeled images

  type(image),allocatable,public :: peeled_image(:)
  ! the peeled images

  type(image),allocatable,dimension(:),public :: peeled_image_sum

  real(dp),allocatable :: d_min(:), d_max(:)
  ! the distance range within which to peeloff photons for each group

  real(dp),allocatable :: column_density(:)

  ! Cache for raytracing

  type cached_source_spectrum
     real(dp),allocatable :: spectrum(:)
  end type cached_source_spectrum

  type(cached_source_spectrum), allocatable, target :: cached_source_spectra(:,:)

  type cached_dust_emissivity
     real(dp),allocatable :: log10_emissivity(:,:)
  end type cached_dust_emissivity

  type(cached_dust_emissivity), allocatable, target :: cached_dust_emissivities(:,:)

  type cached_dust_extinction
     real(dp),allocatable :: chi_nu(:)
  end type cached_dust_extinction

  type(cached_dust_extinction), allocatable, target :: cached_dust_extinctions(:,:)

  type spectrum
     real(dp),allocatable :: spectrum(:)
     real(dp),allocatable :: chi_nu(:)
  end type spectrum

  type(spectrum),allocatable :: tmp(:)

contains

  subroutine peeled_images_adjust_scale(scale)
    implicit none
    real(dp),intent(in) :: scale
    integer :: ig
    do ig=1,n_groups
       call image_scale(peeled_image(ig), scale)
    end do
  end subroutine peeled_images_adjust_scale

  subroutine peeloff_photon(p_orig, polychromatic)

    implicit none

    type(photon),intent(in) :: p_orig
    type(photon) :: p
    real(dp) :: tau
    integer :: ip,ig,iv,id
    type(angle3d_dp) :: a_req, a_sky, a_diff
    type(vector3d_dp) :: v_req
    logical,intent(in) :: polychromatic
    real(dp) :: x_image, y_image
    real(dp) :: tmax
    type(vector3d_dp) :: dr
    real(dp) :: d
    logical :: killed

    do ip=1,n_peeled

       ig = group_id(ip)
       iv = view_id(ip)

       p = p_orig

       ! Set photon properties to pre-interaction (except frequency)
       p%s = p%s_prev
       p%a = p%a_prev
       p%v = p%v_prev

       ! Find peeloff angle
       a_req = a_peeloff(p%r, ip)
       call angle3d_to_vector3d(a_req, v_req)

       ! Find peeloff probability
       if(p%last_isotropic) then
          p%s = stokes_dp(1._dp, 0._dp, 0._dp, 0._dp)
          p%a = a_req
          p%v = v_req
       else
          select case(p%last)
          case('sr')
             call emit_peeloff(p, a_req)
          case('ds','de')
             call interact_peeloff(p, a_req)
          case('me')
             call mrw_peeloff(p, a_req)
          case default
             call error("peeloff_photon","unexpected p%last flag: "//p%last)
          end select
       end if

       ! Update wall and cell that photon is on, because when the direction
       ! of the photon changes, the cell and wall it might be considered on
       ! may change. For example, if a photon is propagaging in the positive x
       ! direction and is in cell 2 on the xmin wall, but gets peeled off in
       ! the negative direction, we need to treat it as if it was in cell 1 on
       ! the xmax wall.
       call place_in_cell(p)

       if(inside_observer(ig)) then
          dr = p%r-r_peeloff(ig)
          d = sqrt(dr.dot.dr)
          tmax = d - d_min(ig)
       else
          d = -(v_req.dot.p%r)
          tmax = - d_min(ig) + d
       end if

       if(d < d_min(ig) .or. d > d_max(ig)) cycle

       if(inside_observer(ig)) then

          if(abs(viewing_angles(ip)%cost) .gt. 1.e-10) then
             call rotate_angle3d(viewing_angles(ip), angle3d_deg(90._dp, 0._dp), a_diff)
             call difference_angle3d(a_diff, -p%a, a_sky)
          else
             a_sky = p%a
          end if

          ! Convert to angles in degrees
          x_image = atan2(a_sky%sinp, a_sky%cosp) * rad2deg
          y_image = atan2(a_sky%sint, a_sky%cost) * rad2deg - 90._dp

          ! Make sure the photon falls inside the image (wrap angles around)
          x_image = peeled_image(ig)%x_max + modulo(x_image - peeled_image(ig)%x_max, 360._dp)
          y_image = peeled_image(ig)%y_min + modulo(y_image - peeled_image(ig)%y_min, 360._dp)

       else

          ! Project onto 2-D plane perpendicular to direction of peeling-off
          dr = p%r - r_peeloff(ig)
          x_image = dr%y * p%a%cosp - dr%x * p%a%sinp
          y_image = dr%z * p%a%sint - dr%y * p%a%cost * p%a%sinp - dr%x * p%a%cost * p%a%cosp

       end if

       if(in_image(peeled_image(ig),x_image, y_image)) then

          if(ignore_optical_depth(ig)) then
             if(polychromatic) then
                column_density = 0._dp
             else
                tau = 0._dp
             end if
             killed = .false.
          else
             if(polychromatic) then
                call grid_escape_column_density(p,tmax,column_density,killed)
             else
                call grid_escape_tau(p,tmax,tau,killed)
             end if
          end if

          ! For inside observer, don't want optical depth to escape grid, just to go to observer!
          ! Need to include 1/d^2!

          if(.not.killed) then

             if(inside_observer(ig)) p%s = p%s / (4._dp * pi * tmax**2._dp)

             if(polychromatic) then

                ! Now the fun begins
                select case(p%emiss_type)
                case(1,2)
                   call get_source_spectrum(p%source_id, ig, tmp(ig)%spectrum)
                case(3)
                   call get_dust_emissivity(p%dust_id, ig, p%emiss_var_id, p%emiss_var_frac, tmp(ig)%spectrum)
                case default
                   stop "unknown emiss_type"
                end select

                tmp(ig)%spectrum = tmp(ig)%spectrum * p%s%i * p%energy

                do id=1,n_dust
                   call get_dust_extinction(id, ig, tmp(ig)%chi_nu)
                   tmp(ig)%spectrum = tmp(ig)%spectrum * exp(- column_density(id) * tmp(ig)%chi_nu)
                end do

                call image_bin_raytraced(peeled_image(ig),p,x_image,y_image,iv, tmp(ig)%spectrum)

             else
                p%s = p%s * exp(-tau)
                call image_bin(peeled_image(ig),p,x_image,y_image,iv)
             end if

          end if

       end if

    end do

  end subroutine peeloff_photon

  subroutine peeled_images_setup(handle, paths, use_raytracing, use_exact_nu, frequencies)

    implicit none

    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: paths(:)

    real(dp),intent(in),optional :: frequencies(:)

    logical,intent(in) :: use_raytracing, use_exact_nu

    integer :: ig,iv,ip = 0

    integer :: n_view

    ! real(dp) :: array_size = 0.

    real(dp),allocatable :: theta(:), phi(:)

    n_groups = size(paths)

    allocate(peeled_image(n_groups))
    allocate(peeled_image_sum(n_groups))
    if(main_process()) write(*,'(" [peeled_images] setting up ",I0," peeled image groups ")') n_groups

    allocate(inside_observer(n_groups))
    allocate(ignore_optical_depth(n_groups))
    allocate(r_peeloff(n_groups))
    allocate(d_min(n_groups))
    allocate(d_max(n_groups))

    ! Find total number of viewing angles
    n_peeled = 0
    do ig=1,n_groups
       call mp_read_keyword(handle, paths(ig), 'inside_observer', inside_observer(ig))
       call mp_read_keyword(handle, paths(ig), 'ignore_optical_depth', ignore_optical_depth(ig))
       call mp_read_keyword(handle, paths(ig), 'n_view', n_view)
       if(.not. n_view > 0) call error("n_view should be a positive integer", "peeled_images_setup")
       call mp_read_keyword(handle, paths(ig), 'd_min', d_min(ig))
       call mp_read_keyword(handle, paths(ig), 'd_max', d_max(ig))
       if(inside_observer(ig).and.d_min(ig) < 0.) then
          if(main_process()) call warn("d_min cannot be < 0. for inside observer - resetting to 0", "peeled_images_setup")
          d_min(ig) = 0.
       end if
       n_peeled = n_peeled + n_view
    end do

    allocate(group_id(n_peeled))
    allocate(view_id(n_peeled))
    allocate(viewing_angles(n_peeled))

    ip = 0

    if(use_raytracing) then
       allocate(cached_source_spectra(n_groups, n_sources))
       allocate(cached_dust_emissivities(n_groups, n_dust))
       allocate(cached_dust_extinctions(n_groups, n_dust))
    end if

    do ig=1,n_groups

       if(inside_observer(ig)) then
          call mp_read_keyword(handle, paths(ig), 'observer_x', r_peeloff(ig)%x)
          call mp_read_keyword(handle, paths(ig), 'observer_y', r_peeloff(ig)%y)
          call mp_read_keyword(handle, paths(ig), 'observer_z', r_peeloff(ig)%z)
       else
          call mp_read_keyword(handle, paths(ig), 'peeloff_x', r_peeloff(ig)%x)
          call mp_read_keyword(handle, paths(ig), 'peeloff_y', r_peeloff(ig)%y)
          call mp_read_keyword(handle, paths(ig), 'peeloff_z', r_peeloff(ig)%z)
       end if

       call mp_read_keyword(handle, paths(ig), 'n_view', n_view)
       call mp_table_read_column_auto(handle, trim(paths(ig))//'/angles', 'theta', theta)
       call mp_table_read_column_auto(handle, trim(paths(ig))//'/angles', 'phi', phi)

       call image_setup(handle,paths(ig),peeled_image(ig),n_view,n_sources,n_dust,use_exact_nu,frequencies)

       if(peeled_image(ig)%use_filters.and.use_raytracing) then
          call error("peeled_images_setup", "filter convolution cannot be used with raytracing")
       end if

       peeled_image(ig)%group_id = ig

       ! If an inside observer, check that the longitudes are inverted
       if(inside_observer(ig)) then
          if(peeled_image(ig)%compute_image) then
             if(peeled_image(ig)%x_min < peeled_image(ig)%x_max) call error("peeled_images_setup", "longitudes should increase towards the left for inside observers")
          end if
          if(peeled_image(ig)%compute_sed) then
            call error("peeled_images_setup", "computing SEDs for inside observers is not supported")
          end if
       end if

       do iv=1,n_view
          ip = ip + 1
          group_id(ip) = ig
          view_id(ip) = iv
          viewing_angles(ip) = angle3d_deg(theta(iv), phi(iv))
       end do

    end do

    allocate(column_density(n_dust))

    allocate(tmp(n_groups))
    do ig=1,n_groups
       allocate(tmp(ig)%spectrum(peeled_image(ig)%n_nu))
       allocate(tmp(ig)%chi_nu(peeled_image(ig)%n_nu))
    end do

  end subroutine peeled_images_setup

  subroutine peeled_images_write(group)

    integer(hid_t),intent(in) :: group
    integer(hid_t) :: g_indiv

    integer :: ig

    character(len=100) :: group_name

    do ig=1,n_groups

       write(group_name, '("group_",I5.5)') ig

       g_indiv = mp_create_group(group, group_name)
       call image_write(peeled_image(ig),g_indiv)
       call mp_write_keyword(g_indiv, '.', 'inside_observer', inside_observer(ig))

       call mp_write_keyword(g_indiv, '.', 'd_min', d_min(ig))
       call mp_write_keyword(g_indiv, '.', 'd_max', d_max(ig))

       call mp_close_group(g_indiv)

    end do

  end subroutine peeled_images_write

  type(angle3d_dp) function a_peeloff(r, ip)
    implicit none
    type(vector3d_dp), intent(in) :: r
    integer, intent(in) :: ip
    integer :: ig
    ig = group_id(ip)
    if(inside_observer(ig)) then
       call vector3d_to_angle3d(r_peeloff(ig)-r, a_peeloff)
    else
       a_peeloff = viewing_angles(ip)
    end if
  end function a_peeloff

  subroutine get_source_spectrum(source_id, group_id, spectrum)

    implicit none

    integer,intent(in) :: source_id, group_id
    real(dp), intent(out) :: spectrum(:)

    ! For now, if the cache is set then we can just use that. In future we
    ! may have to check whether any conditions (e.g. viewing angle) have
    ! changed since last time this was called.

    if(.not.allocated(cached_source_spectra(group_id, source_id)%spectrum)) then

       allocate(cached_source_spectra(group_id,source_id)%spectrum(peeled_image(group_id)%n_nu))

       if(peeled_image(group_id)%use_exact_nu) then
          spectrum = get_spectrum_interp(s(source_id), peeled_image(group_id)%nu)
       else
          spectrum = get_spectrum_binned(s(source_id), peeled_image(group_id)%n_nu, &
               & peeled_image(group_id)%nu_min, peeled_image(group_id)%nu_max)
       end if

       ! Store the cached spectrum
       cached_source_spectra(group_id,source_id)%spectrum = spectrum

    end if

    spectrum = cached_source_spectra(group_id, source_id)%spectrum

  end subroutine get_source_spectrum

  subroutine get_dust_emissivity(dust_id, group_id, jnu_var_id, jnu_var_frac, spectrum)

    implicit none

    integer,intent(in) :: dust_id, group_id, jnu_var_id
    real(dp), intent(in) :: jnu_var_frac
    integer :: ijnu
    real(dp), intent(out) :: spectrum(:)
    real(dp), pointer :: log10_j_nu(:,:)

    ! For now, if the cache is set then we can just use that. In future we
    ! may have to check whether any conditions (e.g. viewing angle) have
    ! changed since last time this was called.

    if(.not.allocated(cached_dust_emissivities(group_id, dust_id)%log10_emissivity)) then

       allocate(cached_dust_emissivities(group_id,dust_id)%log10_emissivity(peeled_image(group_id)%n_nu, size(d(dust_id)%j_nu_var)))

       do ijnu=1,size(d(dust_id)%j_nu_var)

          if(peeled_image(group_id)%use_exact_nu) then
             spectrum = get_j_nu_interp(d(dust_id), peeled_image(group_id)%nu, ijnu)
          else
             spectrum = get_j_nu_binned(d(dust_id), peeled_image(group_id)%n_nu, &
                  & peeled_image(group_id)%nu_min, peeled_image(group_id)%nu_max, ijnu)

          end if

          ! Assign to cache
          cached_dust_emissivities(group_id,dust_id)%log10_emissivity(:,ijnu) = log10(spectrum)

       end do

    end if

    log10_j_nu => cached_dust_emissivities(group_id,dust_id)%log10_emissivity

    spectrum = (log10_j_nu(:,jnu_var_id+1) - log10_j_nu(:,jnu_var_id)) * jnu_var_frac + &
         &      log10_j_nu(:,jnu_var_id)

    spectrum = 10._dp ** spectrum

    where(is_nan(spectrum))
       spectrum = 0.
    end where

  end subroutine get_dust_emissivity

  subroutine get_dust_extinction(dust_id, group_id, extinction)

    implicit none

    integer,intent(in) :: dust_id, group_id
    real(dp), intent(out) :: extinction(:)

    ! For now, if the cache is set then we can just use that.

    if(.not.allocated(cached_dust_extinctions(group_id, dust_id)%chi_nu)) then

       allocate(cached_dust_extinctions(group_id,dust_id)%chi_nu(peeled_image(group_id)%n_nu))

       if(peeled_image(group_id)%use_exact_nu) then
          extinction = get_chi_nu_interp(d(dust_id), peeled_image(group_id)%nu)
       else
          extinction = get_chi_nu_binned(d(dust_id), peeled_image(group_id)%n_nu, &
               & peeled_image(group_id)%nu_min, peeled_image(group_id)%nu_max)

       end if

       ! Store the cached extinction
       cached_dust_extinctions(group_id,dust_id)%chi_nu = extinction

    end if

    extinction = cached_dust_extinctions(group_id, dust_id)%chi_nu

  end subroutine get_dust_extinction

end module peeled_images
