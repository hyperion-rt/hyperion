module peeled_images

  use core_lib
  use mpi_core

  use grid_propagate

  use type_image
  use type_photon

  use sources
  use dust_main
  use dust_interact

  use grid_mrw

  implicit none

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
  type(vector3d_dp), allocatable :: r_peeloff(:)
  type(angle3d_dp), allocatable :: a_peeloff_fixed(:)
  ! angles in which to make peeled images

  type(image),allocatable,public :: peeled_image(:)
  ! the peeled images

  type(image),allocatable,dimension(:),public :: peeled_image_sum

  real(dp),allocatable :: d_min(:), d_max(:)
  ! the distance range within which to peeloff photons for each group

  real(dp),allocatable :: column_density(:)

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
    integer :: ip,ig,iv
    type(angle3d_dp) :: a_req
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

       if(inside_observer(ig)) then
          dr = p%r-r_peeloff(ig)
          d = sqrt(dr.dot.dr)
          tmax = d
       else
          d = -(v_req.dot.p%r)
          tmax = huge(1._dp)
       end if

       if(d < d_min(ig) .or. d > d_max(ig)) return

       if(polychromatic) then
          call grid_escape_column_density(p,tmax,column_density,killed)
       else
          call grid_escape_tau(p,tmax,tau,killed)
       end if

       ! For inside observer, don't want optical depth to escape grid, just to go to observer!
       ! Need to include 1/d^2!

       if(.not.killed) then

          if(inside_observer(ig)) then

             p%s = p%s / (4._dp * pi * tmax**2._dp)

             ! Convert to angles in degrees
             x_image = atan2(p%a%sinp, p%a%cosp) * rad2deg + 180._dp
             y_image = 90._dp - atan2(p%a%sint, p%a%cost) * rad2deg
             if(y_image > 90._dp) then
                y_image = 180._dp - y_image
                x_image = x_image - 180._dp
             end if
             if(x_image < -180._dp) x_image = x_image + 360._dp
             if(x_image > +180._dp) x_image = x_image - 360._dp
          else
             ! Project onto 2-D plane perpendicular to direction of peeling-off
             dr = p%r - r_peeloff(ig)
             x_image = dr%y * p%a%cosp - dr%x * p%a%sinp
             y_image = dr%z * p%a%sint - dr%y * p%a%cost * p%a%sinp - dr%x * p%a%cost * p%a%cosp
          end if

          if(polychromatic) then
             call image_bin_raytraced(peeled_image(ig),p,x_image,y_image,iv,column_density)
          else
             p%s = p%s * exp(-tau)
             call image_bin(peeled_image(ig),p,x_image,y_image,iv)
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

    integer :: ig,id,iv,is,ip = 0

    integer :: n_view

    ! real(dp) :: array_size = 0.

    real(dp),allocatable :: theta(:), phi(:)

    n_groups = size(paths)

    allocate(peeled_image(n_groups))
    allocate(peeled_image_sum(n_groups))
    if(main_process()) write(*,'(" [peeled_images] setting up ",I0," peeled image groups ")') n_groups

    allocate(inside_observer(n_groups))
    allocate(r_peeloff(n_groups))
    allocate(d_min(n_groups))
    allocate(d_max(n_groups))

    ! Find total number of viewing angles
    n_peeled = 0
    do ig=1,n_groups
       call hdf5_read_keyword(handle, paths(ig), 'inside_observer', inside_observer(ig))
       if(inside_observer(ig)) then
          n_view = 1
       else
          call hdf5_read_keyword(handle, paths(ig), 'n_view', n_view)
       end if
       call hdf5_read_keyword(handle, paths(ig), 'd_min', d_min(ig))
       call hdf5_read_keyword(handle, paths(ig), 'd_max', d_max(ig))
       n_peeled = n_peeled + n_view
    end do

    allocate(group_id(n_peeled))
    allocate(view_id(n_peeled))
    allocate(a_peeloff_fixed(n_peeled))

    ip = 0

    do ig=1,n_groups

       if(inside_observer(ig)) then
          n_view = 1
          call hdf5_read_keyword(handle, paths(ig), 'observer_x', r_peeloff(ig)%x)
          call hdf5_read_keyword(handle, paths(ig), 'observer_y', r_peeloff(ig)%y)
          call hdf5_read_keyword(handle, paths(ig), 'observer_z', r_peeloff(ig)%z)
       else
          call hdf5_read_keyword(handle, paths(ig), 'n_view', n_view)
          call hdf5_table_read_column_auto(handle, trim(paths(ig))//'/Angles', 'theta', theta)
          call hdf5_table_read_column_auto(handle, trim(paths(ig))//'/Angles', 'phi', phi)
          call hdf5_read_keyword(handle, paths(ig), 'peeloff_x', r_peeloff(ig)%x)
          call hdf5_read_keyword(handle, paths(ig), 'peeloff_y', r_peeloff(ig)%y)
          call hdf5_read_keyword(handle, paths(ig), 'peeloff_z', r_peeloff(ig)%z)
       end if

       call image_setup(handle,paths(ig),peeled_image(ig),n_view,use_exact_nu,frequencies)

       do iv=1,n_view
          ip = ip + 1
          group_id(ip) = ig
          view_id(ip) = iv
          if(.not.inside_observer(ig)) then
             a_peeloff_fixed(ip) = angle3d_deg(theta(iv), phi(iv))
          else
             a_peeloff_fixed(ip) = angle3d_deg(0._dp,0._dp)
          end if
       end do

       ! Set up raytracing spectra

       if(use_raytracing) then

          call image_raytracing_initialize(peeled_image(ig),n_sources,n_dust)

          do is=1,n_sources
             select case(s(is)%freq_type)
             case(1)
                call image_raytracing_set_spectrum(peeled_image(ig),is,s(is)%spectrum)
             case(2)
                call image_raytracing_set_blackbody(peeled_image(ig),is,s(is)%temperature)
             case default
                print *,"Don't need to set up spectrum for source = ",is
             end select
          end do

          do id=1,n_dust
             call image_raytracing_set_opacity(peeled_image(ig),id,d(id)%nu,d(id)%chi_nu)
             call image_raytracing_set_emissivity(peeled_image(ig),id,d(id)%j_nu_var,d(id)%j_nu)
          end do

       end if

    end do

    allocate(column_density(n_dust))

    ! array_size = array_size + 8 * n_x * n_y * n_w / 1024.**2.   
    ! write(*,'(" [",F6.1,"Mb]")') array_size

  end subroutine peeled_images_setup

  subroutine peeled_images_write(group)

    integer(hid_t),intent(in) :: group
    integer(hid_t) :: g_indiv

    integer :: ig

    character(len=100) :: group_name

    do ig=1,n_groups
       write(group_name, '("Group ",I5.5)') ig
       g_indiv = hdf5_create_group(group, group_name)
       call image_write(peeled_image(ig),g_indiv)
       call hdf5_close_group(g_indiv)
    end do

  end subroutine peeled_images_write

  type(angle3d_dp) function a_peeloff(r, ip)
    implicit none
    type(vector3d_dp), intent(in) :: r
    integer, intent(in) :: ip
    integer :: ig
    ig = group_id(ip)
    if(inside_observer(ig)) then
       call vector3d_to_angle3d(r-r_peeloff(ig), a_peeloff)
    else
       a_peeloff = a_peeloff_fixed(ip)
    end if
  end function a_peeloff

end module peeled_images
