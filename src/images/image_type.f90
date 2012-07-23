! This module is used both for binned and peeloff images

module type_image

  use core_lib
  use mpi_hdf5_io
  use type_photon

  implicit none
  save

  private
  public :: image
  public :: image_scale
  public :: image_bin
  public :: image_bin_raytraced
  public :: image_setup
  public :: image_write
  public :: in_image

  public :: image_raytracing_initialize
  public :: image_raytracing_set_spectrum
  public :: image_raytracing_set_blackbody
  public :: image_raytracing_set_emissivity
  public :: image_raytracing_set_opacity

  public :: source_helper
  public :: dust_helper

  type source_helper
     real(dp), allocatable :: spectrum(:)
  end type source_helper

  type dust_helper
     real(dp), allocatable :: chi(:)
     real(dp), allocatable :: emissivity_var(:)
     real(dp), allocatable :: log10_emissivity(:,:)
  end type dust_helper

  type image

     ! IMAGE PARAMETERS

     integer :: n_x,n_y
     real(dp) :: x_min,x_max,y_min,y_max
     ! number of pixels in each direction, and view limits

     integer :: n_nu
     logical :: use_exact_nu
     integer :: inu_min,inu_max
     real(dp) :: nu_min,nu_max
     real(dp) :: log10_nu_min,log10_nu_max
     real(dp),allocatable :: nu(:)
     ! number of frequencies, and frequency range

     integer :: n_ap
     real(dp) :: ap_min,ap_max
     real(dp) :: log10_ap_min,log10_ap_max
     ! number of apertures for the SED, and aperture range

     integer :: n_view
     type(angle3d_dp), allocatable :: a(:)
     ! number of viewing angles, and viewing angles

     ! IMAGE AND SED

     logical :: compute_image
     type(stokes_dp),allocatable :: img(:,:,:,:,:)
     type(stokes_dp),allocatable :: img2(:,:,:,:,:)
     real(dp),allocatable :: imgn(:,:,:,:,:)
     ! the stokes Image array and the summation array for the variance

     logical :: compute_sed
     type(stokes_dp),allocatable :: sed(:,:,:,:)
     type(stokes_dp),allocatable :: sed2(:,:,:,:)
     real(dp),allocatable :: sedn(:,:,:,:)
     ! the stokes SED array and the summation array for the variance

     ! RAYTRACING

     type(source_helper),allocatable :: sources(:)
     ! spectra of the sources, binned to the image spectral resolution

     type(dust_helper),allocatable :: dust(:)
     ! emissivity of the dust types, binned to the image spectral resolution

     real(dp),allocatable :: tmp_spectrum(:)

     character(len=10) :: track_origin
     logical :: uncertainties

     integer :: io_type

     integer :: n_orig

     integer :: n_sources
     integer :: n_dust

  end type image

  interface image_raytracing_set_spectrum
     module procedure image_raytracing_set_spectrum_pdf
     module procedure image_raytracing_set_spectrum_nu_fnu
  end interface image_raytracing_set_spectrum

contains

  integer function orig(p)
    ! Find the origin flag for a photon
    implicit none
    type(photon),intent(in) :: p
    if(p%scattered) then
       if(p%reprocessed) then
          orig = 4 ! scattered dust emission
       else
          orig = 3 ! scattered source emission
       end if
    else
       if(p%reprocessed) then
          orig = 2 ! dust emission
       else
          orig = 1 ! source emission
       end if
    end if
  end function orig

  subroutine image_scale(img,scale)

    implicit none

    type(image), intent(inout) :: img
    real(dp), intent(in) :: scale

    if(img%compute_image) then
       img%img%i = img%img%i * scale
       img%img%q = img%img%q * scale
       img%img%u = img%img%u * scale
       img%img%v = img%img%v * scale
    end if
    if(img%compute_sed) then
       img%sed%i = img%sed%i * scale
       img%sed%q = img%sed%q * scale
       img%sed%u = img%sed%u * scale
       img%sed%v = img%sed%v * scale
    end if

    if(img%uncertainties) then
       if(img%compute_image) then
          img%img2%i = img%img2%i * scale**2
          img%img2%q = img%img2%q * scale**2
          img%img2%u = img%img2%u * scale**2
          img%img2%v = img%img2%v * scale**2
       end if
       if(img%compute_sed) then
          img%sed2%i = img%sed2%i * scale**2
          img%sed2%q = img%sed2%q * scale**2
          img%sed2%u = img%sed2%u * scale**2
          img%sed2%v = img%sed2%v * scale**2
       end if
    end if

  end subroutine image_scale

  subroutine image_setup(handle,path,img,n_view,n_sources,n_dust,use_exact_nu,frequencies)

    implicit none
    integer(hid_t),intent(in) :: handle
    character(len=*),intent(in) :: path

    type(image),intent(out) :: img
    integer,intent(in) :: n_view, n_sources, n_dust
    logical,intent(in) :: use_exact_nu
    real(dp),intent(in),optional :: frequencies(:)
    real(dp) :: wav_min, wav_max
    integer :: io_bytes

    img%n_view = n_view
    img%n_sources = n_sources
    img%n_dust = n_dust

    call mp_read_keyword(handle, path, 'n_wav',img%n_nu)

    call mp_read_keyword(handle, path, 'compute_image',img%compute_image)
    if(img%compute_image) then
       call mp_read_keyword(handle, path, 'n_x',img%n_x)
       call mp_read_keyword(handle, path, 'n_y',img%n_y)
       call mp_read_keyword(handle, path, 'x_min',img%x_min)
       call mp_read_keyword(handle, path, 'x_max',img%x_max)
       call mp_read_keyword(handle, path, 'y_min',img%y_min)
       call mp_read_keyword(handle, path, 'y_max',img%y_max)
    end if

    call mp_read_keyword(handle, path, 'compute_sed',img%compute_sed)
    if(img%compute_sed) then
       call mp_read_keyword(handle, path, 'n_ap',img%n_ap)
       call mp_read_keyword(handle, path, 'ap_min',img%ap_min)
       call mp_read_keyword(handle, path, 'ap_max',img%ap_max)
       img%log10_ap_min = log10(img%ap_min)
       img%log10_ap_max = log10(img%ap_max)
    end if

    call mp_read_keyword(handle, path, 'track_origin',img%track_origin)

    select case(trim(img%track_origin))
    case('detailed')
       img%n_orig = 2 * (img%n_sources + img%n_dust)
    case('basic', 'yes')
       img%n_orig = 4
    case('no')
       img%n_orig = 1
    case default
       call error("image_setup", "unknown track_origin flag: "//trim(img%track_origin))
    end select

    call mp_read_keyword(handle, path, 'uncertainties',img%uncertainties)

    if(use_exact_nu) then

       img%use_exact_nu = .true.

       call mp_read_keyword(handle, path, 'inu_min',img%inu_min)
       call mp_read_keyword(handle, path, 'inu_max',img%inu_max)

       if(.not.present(frequencies)) call error('image_setup','frequencies should be given if use_exact_nu is .true.')

       allocate(img%nu(size(frequencies)))
       img%nu = frequencies(img%inu_min:img%inu_max)

    else

       img%use_exact_nu = .false.

       call mp_read_keyword(handle, path, 'wav_min',wav_min)
       call mp_read_keyword(handle, path, 'wav_max',wav_max)

       img%nu_min = c_cgs / (wav_max * 1.e-4)
       img%nu_max = c_cgs / (wav_min * 1.e-4)

       img%log10_nu_min = log10(img%nu_min)
       img%log10_nu_max = log10(img%nu_max)

    end if

    ! Allocate arrays

    if(img%compute_image) then
       allocate(img%img(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig))
       if(img%uncertainties) then
          allocate(img%img2(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig))
          allocate(img%imgn(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig))
       end if
    end if

    if(img%compute_sed) then
       allocate(img%sed(img%n_nu,img%n_ap,img%n_view,img%n_orig))
       if(img%uncertainties) then
          allocate(img%sed2(img%n_nu,img%n_ap,img%n_view,img%n_orig))
          allocate(img%sedn(img%n_nu,img%n_ap,img%n_view,img%n_orig))
       end if
    end if

    ! Initialize all arrays

    if(img%compute_image) then
       img%img = stokes_dp(0.,0.,0.,0.)
       if(img%uncertainties) then
          img%img2 = stokes_dp(0.,0.,0.,0.)
          img%imgn = 0._dp
       end if
    end if

    if(img%compute_sed) then
       img%sed = stokes_dp(0.,0.,0.,0.)
       if(img%uncertainties) then
          img%sed2 = stokes_dp(0.,0.,0.,0.)
          img%sedn = 0._dp
       end if
    end if

    call mp_read_keyword(handle, path, 'io_bytes',io_bytes)

    select case(io_bytes)
    case(4)
       img%io_type = sp
    case(8)
       img%io_type = dp
    case default
       call error("image_setup", "unexpected value of io_bytes (should be 4 or 8)")
    end select

  end subroutine image_setup

  subroutine find_sed_bin(img,x_image,y_image,ir)

    implicit none

    type(image),intent(in) :: img  ! Image
    real(dp),intent(in) :: x_image, y_image
    integer,intent(out) :: ir
    real(dp) :: log10_r_image

    log10_r_image = log10(sqrt(x_image*x_image + y_image*y_image))

    if(log10_r_image < img%log10_ap_min .or. img%n_ap == 1) then
       ir = 1
    else
       ir = ipos(img%log10_ap_min,img%log10_ap_max,log10_r_image,img%n_ap-1) + 1
    end if

  end subroutine find_sed_bin

  subroutine find_image_bin(img,x_image,y_image,ix,iy)

    implicit none

    type(image),intent(in)    :: img  ! Image
    real(dp),intent(in) :: x_image, y_image
    integer,intent(out) :: ix,iy

    ix = ipos(img%x_min,img%x_max,x_image,img%n_x)
    iy = ipos(img%y_min,img%y_max,y_image,img%n_y)

  end subroutine find_image_bin

  logical function in_image(img, x_image, y_image, nu)

    implicit none

    type(image),intent(in) :: img
    real(dp),intent(in) :: x_image, y_image
    real(dp),intent(in),optional :: nu

    in_image = .false.

    ! Check if photon would end up inside the frequency range. Note that we
    ! can only know for sure that the photon is outside the range in which
    ! case in_image is already set to .false. so we just return. If it's
    ! inside the frequency range, then we still have to check if it would
    ! fall in the image or SEDs.
    if(present(nu)) then
       if(nu < img%nu_min .or. nu > img%nu_max) return
    end if

    ! Check if photon would fall in image
    if(img%compute_image) then
       if((x_image >= img%x_min .and. x_image <= img%x_max) .or. (x_image <= img%x_min .and. x_image >= img%x_max)) then
          if((y_image >= img%y_min .and. y_image <= img%y_max) .or. (y_image <= img%y_min .and. y_image >= img%y_max)) then
             in_image = .true.
             return
          end if
       end if
    end if

    ! Check if photon would fall in SED
    if(img%compute_sed) then
       if(x_image*x_image + y_image*y_image <= img%ap_max * img%ap_max) then
          in_image = .true.
          return
       end if
    end if

  end function in_image

  subroutine image_bin(img,p,x_image,y_image,im)

    implicit none

    type(image),intent(inout)    :: img  ! Image
    type(photon),intent(in) :: p
    real(dp),intent(in) :: x_image, y_image
    integer,intent(in) :: im ! sub-image to bin into
    real(dp) :: log10_nu_image
    integer :: ix,iy,ir,inu,io ! Bins

    if(img%compute_image.and..not.allocated(img%img)) call error('bin_photon','Image not allocated')
    if(img%compute_sed.and..not.allocated(img%sed)) call error('bin_photon','SED not allocated')

    if(is_nan(p%energy)) then
       call warn("image_bin","photon has NaN energy - ignoring")
       return
    end if

    if(is_nan(p%s%i)) then
       call warn("image_bin","photon has NaN flux - ignoring")
       return
    end if

    ! Find frequency bin
    if(img%use_exact_nu) then
       inu = p%inu - img%inu_min + 1
    else
       log10_nu_image = log10(p%nu)
       inu = ipos(img%log10_nu_min,img%log10_nu_max,log10_nu_image,img%n_nu)
    end if

    ! Find origin flag
    if(trim(img%track_origin) == 'detailed') then
       io = ((orig(p) - mod(orig(p),2)) * img%n_sources + (orig(p) - mod(orig(p)+1,2) - 1) * img%n_dust) / 2
       if(mod(orig(p),2)==0) then
          io = io + p%dust_id
       else
          io = io + p%source_id
       end if
    else if(trim(img%track_origin) == 'basic') then
       io = orig(p)
    else
       io = 1
    end if

    if(inu >= 1 .and. inu <= img%n_nu) then
       if(img%compute_image) then
          call find_image_bin(img,x_image,y_image,ix,iy)
          if(ix >= 1 .and. ix <= img%n_x) then
             if(iy >= 1 .and. iy <= img%n_y) then
                img%img(inu,ix,iy,im,io) = img%img(inu,ix,iy,im,io) + p%s * p%energy
                if(img%uncertainties) then
                   img%img2(inu,ix,iy,im,io) = img%img2(inu,ix,iy,im,io) + p%s**2._dp * p%energy**2._dp
                   img%imgn(inu,ix,iy,im,io) = img%imgn(inu,ix,iy,im,io) + 1._dp
                end if
             end if
          end if
       end if
       if(img%compute_sed) then
          call find_sed_bin(img,x_image,y_image,ir)
          if(ir >= 1 .and. ir <= img%n_ap) then
             img%sed(inu,ir,im,io) = img%sed(inu,ir,im,io) + p%s * p%energy
             if(img%uncertainties) then
                img%sed2(inu,ir,im,io) = img%sed2(inu,ir,im,io) + p%s**2._dp * p%energy**2._dp
                img%sedn(inu,ir,im,io) = img%sedn(inu,ir,im,io) + 1._dp
             end if
          end if
       end if
    end if

  end subroutine image_bin

  subroutine image_bin_raytraced(img,p,x_image,y_image,im,column_density)

    implicit none

    type(image),intent(inout)    :: img  ! Image
    type(photon),intent(in) :: p
    real(dp),intent(in) :: x_image, y_image
    integer,intent(in) :: im ! sub-image to bin into
    real(dp),intent(in) :: column_density(:)
    integer :: ix,iy,ir,iw,id,io ! Bins

    if(img%compute_image.and..not.allocated(img%img)) call error('bin_photon','Image not allocated')
    if(img%compute_sed.and..not.allocated(img%sed)) call error('bin_photon','SED not allocated')

    if(is_nan(p%energy)) then
       call warn("image_bin_raytraced","photon has NaN energy - ignoring")
       return
    end if

    if(is_nan(p%s%i)) then
       call warn("image_bin_raytraced","photon has NaN flux - ignoring")
       return
    end if

    ! Find origin flag
    if(trim(img%track_origin) == 'detailed') then
       io = ((orig(p) - mod(orig(p),2)) * img%n_sources + (orig(p) - mod(orig(p)+1,2) - 1) * img%n_dust) / 2
       if(mod(orig(p),2)==0) then
          io = io + p%dust_id
       else
          io = io + p%source_id
       end if
    else if(trim(img%track_origin) == 'basic') then
       io = orig(p)
    else
       io = 1
    end if

    ! Now the fun begins
    select case(p%emiss_type)
    case(1,2)
       img%tmp_spectrum = img%sources(p%source_id)%spectrum
    case(3)
       img%tmp_spectrum = (img%dust(p%dust_id)%log10_emissivity(:,p%emiss_var_id+1) - &
            &                   img%dust(p%dust_id)%log10_emissivity(:,p%emiss_var_id)) * p%emiss_var_frac + &
            &                   img%dust(p%dust_id)%log10_emissivity(:,p%emiss_var_id)
       img%tmp_spectrum = 10._dp**(img%tmp_spectrum)
       where(is_nan(img%tmp_spectrum))
          img%tmp_spectrum = 0.
       end where
    case default
       stop "unknown emiss_type"
    end select

    img%tmp_spectrum = img%tmp_spectrum * p%s%i * p%energy

    do id=1,size(column_density)
       img%tmp_spectrum = img%tmp_spectrum * exp(- column_density(id) * img%dust(id)%chi)
    end do

    if(img%compute_image) then
       call find_image_bin(img,x_image,y_image,ix,iy)
       if(ix >= 1 .and. ix <= img%n_x) then
          if(iy >= 1 .and. iy <= img%n_y) then
             do iw=1,img%n_nu
                img%img(iw,ix,iy,im,io)%i = img%img(iw,ix,iy,im,io)%i + img%tmp_spectrum(iw)
                if(img%uncertainties) then
                   img%img2(iw,ix,iy,im,io)%i = img%img2(iw,ix,iy,im,io)%i + img%tmp_spectrum(iw)**2._dp
                   img%imgn(iw,ix,iy,im,io) = img%imgn(iw,ix,iy,im,io) + 1._dp
                end if
             end do
          end if
       end if
    end if

    if(img%compute_sed) then
       call find_sed_bin(img,x_image,y_image,ir)
       if(ir >= 1 .and. ir <= img%n_ap) then
          do iw=1,img%n_nu
             img%sed(iw,ir,im,io)%i = img%sed(iw,ir,im,io)%i + img%tmp_spectrum(iw)
             if(img%uncertainties) then
                img%sed2(iw,ir,im,io)%i = img%sed2(iw,ir,im,io)%i + img%tmp_spectrum(iw)**2._dp
                img%sedn(iw,ir,im,io) = img%sedn(iw,ir,im,io) + 1._dp
             end if
          end do
       end if
    end if

  end subroutine image_bin_raytraced

  subroutine image_write(img,group)

    implicit none

    type(image),intent(in)      :: img  ! Image to write out
    integer(hid_t),intent(in) :: group

    real(dp),allocatable :: cube5d(:,:,:,:,:)
    real(dp),allocatable :: cube5de(:,:,:,:,:)
    real(dp),allocatable :: cube6d(:,:,:,:,:,:)
    real(dp),allocatable :: cube6de(:,:,:,:,:,:)

    integer  :: ia,inu

    real(dp) :: dnunorm

    ! The factor for SED normalization is derived in the following way:
    !
    ! Given n_nu equal bins in log10 space from nu_min to nu_max, the center
    ! of the bins is given by:
    !
    ! log10[nu(j)] = [log10(nu_max) - log10(nu_min)] * (j - 0.5) / n_nu
    !
    ! Thus:
    !
    ! nu(j) = (nu_max / nu_min) ** ((j - 0.5) / n_nu)
    !
    ! Therefore, the width of bins is given by
    !
    ! dnu(j) = nu(j + 0.5) - nu(j - 0.5)
    ! dnu(j) = (nu_max / nu_min) ** (j / n_nu)
    !        - (nu_max / nu_min) ** ((j - 1) / n_nu)
    !
    ! The SED normalization is done by converting F_nu * dnu to nu * Fnu, i.e.
    ! dividing by dnu(j) and multiplying by nu(j), i.e. dividing by dnunorm
    ! where:
    !
    ! dnunorm = dnu(j) / nu(j)
    ! dnunorm = (nu_max / nu_min) ** (+0.5 / n_nu)
    !         - (nu_max / nu_min) ** (-0.5 / n_nu)

    dnunorm = (img%nu_max / img%nu_min) ** (+0.5_dp / real(img%n_nu, dp)) &
         - (img%nu_max / img%nu_min) ** (-0.5_dp / real(img%n_nu, dp))

    if(img%compute_sed) then

       write(*,'(" [image_write] writing out SEDs")')

       allocate(cube5d(img%n_nu, img%n_ap, img%n_view,img%n_orig,4))
       if(img%uncertainties) allocate(cube5de(img%n_nu, img%n_ap, img%n_view,img%n_orig,4))

       cube5d(:,:,:,:,1) = img%sed%i
       cube5d(:,:,:,:,2) = img%sed%q
       cube5d(:,:,:,:,3) = img%sed%u
       cube5d(:,:,:,:,4) = img%sed%v

       if(img%uncertainties) then
          where(img%sedn > 1)
             cube5de(:,:,:,:,1) = sqrt((img%sed2%i + (img%sed%i)**2 / img%sedn) / (img%sedn - 1)) * sqrt(img%sedn)
             cube5de(:,:,:,:,2) = sqrt((img%sed2%q + (img%sed%q)**2 / img%sedn) / (img%sedn - 1)) * sqrt(img%sedn)
             cube5de(:,:,:,:,3) = sqrt((img%sed2%u + (img%sed%u)**2 / img%sedn) / (img%sedn - 1)) * sqrt(img%sedn)
             cube5de(:,:,:,:,4) = sqrt((img%sed2%v + (img%sed%v)**2 / img%sedn) / (img%sedn - 1)) * sqrt(img%sedn)
          elsewhere
             cube5de(:,:,:,:,1) = 0._dp
             cube5de(:,:,:,:,2) = 0._dp
             cube5de(:,:,:,:,3) = 0._dp
             cube5de(:,:,:,:,4) = 0._dp
          end where
       end if

       if(.not.img%use_exact_nu) then
          cube5d = cube5d / dnunorm
          if(img%uncertainties) cube5de = cube5de / dnunorm
       else
          do inu=1,img%n_nu
             cube5d(inu,:,:,:,:) = cube5d(inu,:,:,:,:) * img%nu(inu)
             if(img%uncertainties) cube5de(inu,:,:,:,:) = cube5de(inu,:,:,:,:) * img%nu(inu)
          end do
       end if

       do ia=2,img%n_ap
          cube5d(:,ia,:,:,:) = cube5d(:,ia-1,:,:,:) + cube5d(:,ia,:,:,:)
          if(img%uncertainties) cube5de(:,ia,:,:,:) = sqrt(cube5de(:,ia-1,:,:,:)**2 + cube5de(:,ia,:,:,:)**2)
       end do

       select case(img%io_type)
       case(sp)
          call mp_write_array(group, 'seds', real(cube5d, sp))
          if(img%uncertainties) call mp_write_array(group, 'seds_unc', real(cube5de, sp))
       case(dp)
          call mp_write_array(group, 'seds', real(cube5d, dp))
          if(img%uncertainties) call mp_write_array(group, 'seds_unc', real(cube5de, dp))
       case default
          call error("image_write","unexpected value of img%io_type (should be sp or dp)")
       end select

       if(.not.img%use_exact_nu) then
          call mp_write_keyword(group, 'seds','numin',img%nu_min)
          call mp_write_keyword(group, 'seds','numax',img%nu_max)
       end if
       call mp_write_keyword(group, 'seds','apmin',img%ap_min)
       call mp_write_keyword(group, 'seds','apmax',img%ap_max)

       call mp_write_keyword(group, 'seds', 'track_origin', img%track_origin)
       if(trim(img%track_origin) == 'detailed') then
          call mp_write_keyword(group, 'seds', 'n_sources', img%n_sources)
          call mp_write_keyword(group, 'seds', 'n_dust', img%n_dust)
       end if

    end if

    if(img%compute_image) then

       write(*,'(" [image_write] writing out images")')

       allocate(cube6d(img%n_nu, img%n_x, img%n_y, img%n_view, img%n_orig, 4))
       if(img%uncertainties) allocate(cube6de(img%n_nu, img%n_x, img%n_y, img%n_view, img%n_orig, 4))

       cube6d(:,:,:,:,:,1) = img%img%i
       cube6d(:,:,:,:,:,2) = img%img%q
       cube6d(:,:,:,:,:,3) = img%img%u
       cube6d(:,:,:,:,:,4) = img%img%v

       if(img%uncertainties) then
          where(img%imgn > 1)
             cube6de(:,:,:,:,:,1) = sqrt((img%img2%i + (img%img%i)**2 / img%imgn) / (img%imgn - 1)) * sqrt(img%imgn)
             cube6de(:,:,:,:,:,2) = sqrt((img%img2%q + (img%img%q)**2 / img%imgn) / (img%imgn - 1)) * sqrt(img%imgn)
             cube6de(:,:,:,:,:,3) = sqrt((img%img2%u + (img%img%u)**2 / img%imgn) / (img%imgn - 1)) * sqrt(img%imgn)
             cube6de(:,:,:,:,:,4) = sqrt((img%img2%v + (img%img%v)**2 / img%imgn) / (img%imgn - 1)) * sqrt(img%imgn)
          elsewhere
             cube6de(:,:,:,:,:,1) = 0._dp
             cube6de(:,:,:,:,:,2) = 0._dp
             cube6de(:,:,:,:,:,3) = 0._dp
             cube6de(:,:,:,:,:,4) = 0._dp
          end where
       end if

       if(.not.img%use_exact_nu) then
          cube6d = cube6d / dnunorm
          if(img%uncertainties) cube6de = cube6de / dnunorm
       else
          do inu=1,img%n_nu
             cube6d(inu,:,:,:,:,:) = cube6d(inu,:,:,:,:,:) * img%nu(inu)
             if(img%uncertainties) cube6de(inu,:,:,:,:,:) = cube6de(inu,:,:,:,:,:) * img%nu(inu)
          end do
       end if

       select case(img%io_type)
       case(sp)
          call mp_write_array(group, 'images', real(cube6d, sp))
          if(img%uncertainties) call mp_write_array(group, 'images_unc', real(cube6de, sp))
       case(dp)
          call mp_write_array(group, 'images', real(cube6d, dp))
          if(img%uncertainties) call mp_write_array(group, 'images_unc', real(cube6de, dp))
       case default
          call error("image_write","unexpected value of img%io_type (should be sp or dp)")
       end select

       if(.not.img%use_exact_nu) then
          call mp_write_keyword(group, 'images','numin',img%nu_min)
          call mp_write_keyword(group, 'images','numax',img%nu_max)
       end if
       call mp_write_keyword(group, 'images','xmin',img%x_min)
       call mp_write_keyword(group, 'images','xmax',img%x_max)
       call mp_write_keyword(group, 'images','ymin',img%y_min)
       call mp_write_keyword(group, 'images','ymax',img%y_max)

       call mp_write_keyword(group, 'images', 'track_origin', img%track_origin)
       if(trim(img%track_origin) == 'detailed') then
          call mp_write_keyword(group, 'images', 'n_sources', img%n_sources)
          call mp_write_keyword(group, 'images', 'n_dust', img%n_dust)
       end if

    end if

    if(img%use_exact_nu) then
       call mp_table_write_header(group, 'frequencies',img%inu_max - img%inu_min + 1,1,(/'nu'/),(/1/),(/h5t_ieee_f64le/))
       call mp_table_write_column(group, 'frequencies','nu',img%nu(img%inu_min:img%inu_max))
    end if

    write(*,'(" [image_write] done")')

  end subroutine image_write

  ! RAYTRACING

  subroutine image_raytracing_initialize(img,n_sources,n_dust)
    implicit none
    type(image), intent(inout) :: img
    integer,intent(in) :: n_sources, n_dust
    allocate(img%sources(n_sources))
    allocate(img%dust(n_dust))
    allocate(img%tmp_spectrum(img%n_nu))
  end subroutine image_raytracing_initialize

  subroutine image_raytracing_set_spectrum_pdf(img,source_id,spectrum)

    implicit none

    type(image), intent(inout) :: img
    integer,intent(in) :: source_id
    type(pdf_dp),intent(in) :: spectrum

    call image_raytracing_set_spectrum(img,source_id,spectrum%x,spectrum%pdf)

  end subroutine image_raytracing_set_spectrum_pdf

  subroutine image_raytracing_set_spectrum_nu_fnu(img,source_id,nu,fnu)

    implicit none

    type(image), intent(inout) :: img
    integer,intent(in) :: source_id
    real(dp),intent(in) :: nu(:), fnu(:)
    real(dp), dimension(img%n_nu) :: nu_new, e_new
    real(dp) :: numin, numax
    integer :: inu


    allocate(img%sources(source_id)%spectrum(img%n_nu))

    if(img%use_exact_nu) then

       img%sources(source_id)%spectrum = interp1d_loglog(nu, fnu, img%nu, bounds_error=.false., fill_value=0._dp)

    else

       do inu=1,img%n_nu

          numin = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu-1, dp) / real(img%n_nu, dp))
          numax = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu, dp) / real(img%n_nu, dp))

          nu_new(inu) = 10._dp**(0.5_dp * (log10(numin) + log10(numax)))
          e_new(inu) = integral_loglog(nu, fnu, numin, numax)

       end do

       e_new = e_new / integral_loglog(nu, fnu)

       img%sources(source_id)%spectrum = e_new

    end if

  end subroutine image_raytracing_set_spectrum_nu_fnu

  subroutine image_raytracing_set_blackbody(img,source_id,temperature)

    implicit none

    type(image), intent(inout) :: img
    integer,intent(in) :: source_id
    real(dp),intent(in) :: temperature
    integer :: inu

    integer :: n_nu
    real(dp),allocatable :: nu(:)

    real(dp) :: log10_nu_min, log10_nu_max

    if(img%use_exact_nu) then

       call image_raytracing_set_spectrum(img,source_id,img%nu,normalized_B_nu(img%nu,temperature))

    else

       log10_nu_min = log10(3.e9_dp)  ! 10 cm (0.05K)
       log10_nu_max = log10(3.e16_dp) ! 10 nm (1/2 million K)

       n_nu = ceiling((log10_nu_max - log10_nu_min) * 100000)
       allocate(nu(n_nu))
       do inu=1,n_nu
          nu(inu) = 10._dp**(real(inu-1,dp)/real(n_nu-1,dp)*(log10_nu_max - log10_nu_min) + log10_nu_min)
       end do

       call image_raytracing_set_spectrum(img,source_id,nu,normalized_B_nu(nu,temperature))

    end if

  contains

    elemental real(dp) function normalized_B_nu(nu,T)
      implicit none
      real(dp),intent(in) :: nu,T
      real(dp),parameter :: a = two * h_cgs / c_cgs / c_cgs / stef_boltz * pi
      real(dp),parameter :: b = h_cgs / k_cgs
      real(dp) :: T4
      T4 = T*T*T*T
      normalized_B_nu = a * nu * nu * nu / ( exp(b*nu/T) - one) / T4
    end function normalized_B_nu

  end subroutine image_raytracing_set_blackbody

  subroutine image_raytracing_set_emissivity(img,dust_id,emissivity_var,emissivity)

    implicit none

    type(image),intent(inout) :: img
    integer,intent(in) :: dust_id
    real(dp),intent(in) :: emissivity_var(:)
    type(pdf_dp),intent(in) :: emissivity(:)
    real(dp), dimension(img%n_nu) :: nu_new
    real(dp), dimension(img%n_nu, size(emissivity_var)) :: e_new
    real(dp) :: numin, numax
    integer :: inu
    integer :: j

    if(img%use_exact_nu) then

       do j=1,size(emissivity_var)
          e_new(:,j) = interp1d_loglog(emissivity(j)%x, emissivity(j)%pdf, img%nu, bounds_error=.false., fill_value=0._dp)
       end do

    else

       do j=1,size(emissivity_var)

          do inu=1,img%n_nu

             numin = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu-1, dp) / real(img%n_nu, dp))
             numax = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu, dp) / real(img%n_nu, dp))

             nu_new(inu) = 10._dp**(0.5_dp * (log10(numin) + log10(numax)))
             e_new(inu,j) = integral_loglog(emissivity(j)%x, emissivity(j)%pdf, numin, numax)

          end do

          e_new(:,j) = e_new(:,j) / integral_loglog(emissivity(j)%x, emissivity(j)%pdf)

       end do

    end if

    allocate(img%dust(dust_id)%emissivity_var(size(emissivity_var)))
    img%dust(dust_id)%emissivity_var = emissivity_var

    allocate(img%dust(dust_id)%log10_emissivity(img%n_nu,size(emissivity_var)))
    do j=1,size(emissivity_var)
       img%dust(dust_id)%log10_emissivity(:,j) = log10(e_new(:,j))
    end do

  end subroutine image_raytracing_set_emissivity

  subroutine image_raytracing_set_opacity(img,dust_id,nu,chi_nu)

    implicit none

    type(image),intent(inout) :: img
    integer,intent(in) :: dust_id
    real(dp),intent(in) :: nu(:), chi_nu(:)
    real(dp), dimension(img%n_nu) :: nu_new, chi_nu_new
    real(dp) :: numin, numax
    integer :: inu

    if(img%use_exact_nu) then

       chi_nu_new = interp1d_loglog(nu, chi_nu, img%nu, bounds_error=.false., fill_value=0._dp)

    else

       do inu=1,img%n_nu

          numin = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu-1, dp) / real(img%n_nu, dp))
          numax = 10._dp**(img%log10_nu_min + (img%log10_nu_max - img%log10_nu_min) * real(inu, dp) / real(img%n_nu, dp))

          nu_new(inu) = 10._dp**(0.5_dp * (log10(numin) + log10(numax)))
          chi_nu_new(inu) = integral_loglog(nu, chi_nu, numin, numax) / (numax - numin)

       end do

    end if

    allocate(img%dust(dust_id)%chi(img%n_nu))
    img%dust(dust_id)%chi = chi_nu_new

  end subroutine image_raytracing_set_opacity

end module type_image
