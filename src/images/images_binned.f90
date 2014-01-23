module binned_images

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use type_image
  use type_photon

  use sources, only : n_sources
  use dust_main, only : n_dust

  implicit none
  save

  private

  public :: binned_images_setup
  public :: binned_images_bin_photon
  public :: binned_images_write
  public :: binned_images_adjust_scale

  logical,public :: make_binned_images

  integer,public :: n_phi,n_theta
  ! number of bins in the phi and theta direction

  type(image),public :: binned_image
  ! the binned images

  type(image),public :: binned_image_sum

contains

  subroutine binned_images_adjust_scale(scale)
    implicit none
    real(dp),intent(in) :: scale
    call image_scale(binned_image, scale*real(n_theta, dp)*real(n_phi, dp))
  end subroutine binned_images_adjust_scale

  subroutine binned_images_setup(handle, path)

    implicit none

    integer(hid_t) :: handle
    character(len=*),intent(in) :: path

    call mp_read_keyword(handle, path, 'n_phi', n_phi)
    call mp_read_keyword(handle, path, 'n_theta', n_theta)

    if(main_process()) write(*,'(" [binned_images] setting up ",I0," binned images ")') n_theta*n_phi

    call image_setup(handle, path ,binned_image,n_theta*n_phi,n_sources,n_dust,.false.)

  end subroutine binned_images_setup

  subroutine binned_images_bin_photon(p)

    implicit none

    type(photon),intent(in) :: p

    real(dp) :: phi

    integer :: it,ip
    real(dp) :: x_image, y_image


    phi = atan2(p%a%sinp,p%a%cosp)
    if(phi < 0.) phi = phi + twopi

    it = ipos(-1._dp,+1._dp,p%a%cost,n_theta)
    ip = ipos(0._dp,twopi,phi,n_phi)

    ! Project onto 2-D plane perpendicular to direction of travel
    x_image = p%r%y * p%a%cosp - p%r%x * p%a%sinp
    y_image = p%r%z * p%a%sint - p%r%y * p%a%cost * p%a%sinp - p%r%x * p%a%cost * p%a%cosp

    call image_bin(binned_image,p,x_image,y_image,image_id(it,ip))

  end subroutine binned_images_bin_photon

  integer function image_id(it, ip)
    implicit none
    integer,intent(in) :: it, ip
    image_id = n_phi * (it-1) + ip
  end function image_id

  subroutine binned_images_write(group)
    implicit none
    integer(hid_t),intent(in) :: group
    call image_write(binned_image,group)
  end subroutine binned_images_write

end module binned_images
