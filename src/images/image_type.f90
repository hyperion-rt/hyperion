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

  type filter

     real(dp) :: n_wav
     real(dp), allocatable :: nu(:), tr(:)
     real(dp) :: nu0

  end type filter

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

     ! Here we used to use the type `stokes_dp` type for the image/SED arrays,
     ! but this turns out to not be ideal because there is no way to turn off
     ! the stokes component. However, if instead we simply use the last
     ! dimension of the array to represent the stokes vector, we can allocate
     ! the array differently depending on whether we care about the
     ! polarization elements.

     integer :: n_stokes = 4

     logical :: compute_image
     real(dp),allocatable :: img(:,:,:,:,:,:)
     real(dp),allocatable :: img2(:,:,:,:,:,:)
     real(dp),allocatable :: imgn(:,:,:,:,:,:)
     ! the stokes Image array and the summation array for the variance

     logical :: compute_sed
     real(dp),allocatable :: sed(:,:,:,:,:)
     real(dp),allocatable :: sed2(:,:,:,:,:)
     real(dp),allocatable :: sedn(:,:,:,:,:)
     ! the stokes SED array and the summation array for the variance

     ! RAYTRACING

     type(dust_helper),allocatable :: dust(:)
     ! emissivity of the dust types, binned to the image spectral resolution

     real(dp),allocatable :: tmp_stokes(:)

     character(len=11) :: track_origin
     integer :: track_n_scat
     logical :: uncertainties

     integer :: io_type

     integer :: n_orig

     integer :: n_sources
     integer :: n_dust

     integer :: group_id = 999

     logical :: use_filters
     type(filter), allocatable :: filters(:)

  end type image

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

    if(img%compute_image) img%img = img%img * scale
    if(img%compute_sed) img%sed = img%sed * scale

    if(img%uncertainties) then
       if(img%compute_image) img%img2 = img%img2 * scale**2
       if(img%compute_sed) img%sed2 = img%sed2 * scale**2
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
    logical :: compute_stokes
    integer :: ig, n_filt
    character(len=100) :: group_name

    img%n_view = n_view
    img%n_sources = n_sources
    img%n_dust = n_dust

    ! Set up filters
    if(mp_exists_keyword(handle, path, 'use_filters')) then
       call mp_read_keyword(handle, path, 'use_filters',img%use_filters)
       if(img%use_filters) then
          if(use_exact_nu) call error("image_setup", "cannot use filters in monochromatic mode")
       end if
       call mp_read_keyword(handle, path, 'n_filt',img%n_nu)
    else
       img%use_filters = .false.
       call mp_read_keyword(handle, path, 'n_wav',img%n_nu)
    end if

    if(img%n_nu < 1) call error("image_setup", "n_nu should be >= 1")

    if(mp_exists_keyword(handle, path, 'compute_stokes')) then
       call mp_read_keyword(handle, path, 'compute_stokes',compute_stokes)
    else
       compute_stokes = .true.
    end if

    if(compute_stokes) then
       img%n_stokes = 4
    else
       img%n_stokes = 1
    end if

    allocate(img%tmp_stokes(img%n_stokes))

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

    if(mp_exists_keyword(handle, path, 'track_n_scat')) then
       call mp_read_keyword(handle, path, 'track_n_scat',img%track_n_scat)
    else
       img%track_n_scat = 0
    end if

    select case(trim(img%track_origin))
    case('scatterings')
       img%n_orig = 4 + 2 * img%track_n_scat
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
       if(img%inu_min < 1 .or. img%inu_min > size(frequencies)) call error('image_setup','inu_min value is out of range')
       if(img%inu_max < 1 .or. img%inu_max > size(frequencies)) call error('image_setup','inu_max value is out of range')

       allocate(img%nu(img%n_nu))

       img%nu = frequencies(img%inu_min:img%inu_max)

       if(img%n_nu /= size(img%nu)) call error('image_setup','n_nu should match length of frequencies array')

    else

       img%use_exact_nu = .false.

       if(.not.img%use_filters) then

          call mp_read_keyword(handle, path, 'wav_min',wav_min)
          call mp_read_keyword(handle, path, 'wav_max',wav_max)

          img%nu_min = c_cgs / (wav_max * 1.e-4)
          img%nu_max = c_cgs / (wav_min * 1.e-4)

          img%log10_nu_min = log10(img%nu_min)
          img%log10_nu_max = log10(img%nu_max)

       else

          allocate(img%filters(img%n_nu))
          do ig=1,img%n_nu
             write(group_name, '("filter_",I5.5)') ig
             call mp_table_read_column_auto(handle, trim(path)//"/"//group_name, 'nu', img%filters(ig)%nu)
             call mp_table_read_column_auto(handle, trim(path)//"/"//group_name, 'tn', img%filters(ig)%tr)
             call mp_read_keyword(handle, trim(path)//"/"//group_name, 'nu0', img%filters(ig)%nu0)
          end do
       end if

    end if

    ! Allocate arrays

    if(img%compute_image) then
       allocate(img%img(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig,img%n_stokes))
       if(img%uncertainties) then
          allocate(img%img2(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig,img%n_stokes))
          allocate(img%imgn(img%n_nu,img%n_x,img%n_y,img%n_view,img%n_orig,img%n_stokes))
       end if
    end if

    if(img%compute_sed) then
       allocate(img%sed(img%n_nu,img%n_ap,img%n_view,img%n_orig,img%n_stokes))
       if(img%uncertainties) then
          allocate(img%sed2(img%n_nu,img%n_ap,img%n_view,img%n_orig,img%n_stokes))
          allocate(img%sedn(img%n_nu,img%n_ap,img%n_view,img%n_orig,img%n_stokes))
       end if
    end if

    ! Initialize all arrays

    if(img%compute_image) then
       img%img = 0._dp
       if(img%uncertainties) then
          img%img2 = 0._dp
          img%imgn = 0._dp
       end if
    end if

    if(img%compute_sed) then
       img%sed = 0._dp
       if(img%uncertainties) then
          img%sed2 = 0._dp
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
    integer :: inu,io ! Bins
    integer :: iorig, ifilt
    real(dp) :: transmission

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
       iorig = orig(p)
       io = ((iorig - mod(iorig,2)) * img%n_sources + (iorig - mod(iorig+1,2) - 1) * img%n_dust) / 2
       if(mod(iorig,2)==0) then
          io = io + p%dust_id
       else
          io = io + p%source_id
       end if
    else if(trim(img%track_origin) == 'scatterings') then
       ! The first slice stores photons with no scattering, then slice 2, 3,
       ! etc. store the first scatterings, then a slice is added for the
       ! remaining flux.
       if(p%n_scat > img%track_n_scat) then
          io = img%track_n_scat + 2
       else
          io = p%n_scat + 1
       end if
       if(p%reprocessed) io = io + (img%track_n_scat + 2)
    else if(trim(img%track_origin) == 'basic') then
       io = orig(p)
    else
       io = 1
    end if

    if(img%use_filters) then
       do ifilt=1,size(img%filters)
          transmission = interp1d(img%filters(ifilt)%nu,&
               &                  img%filters(ifilt)%tr,&
               &                  p%nu,bounds_error=.false., fill_value=0._dp)
          if(transmission > 0._dp) then
             call image_bin_single(img, p, x_image, y_image, im, ifilt, io, transmission)
          end if
       end do
    else
       call image_bin_single(img, p, x_image, y_image, im, inu, io, 1._dp)
    end if

  end subroutine image_bin

  subroutine image_bin_single(img,p,x_image,y_image,im,inu,io, transmission)

    implicit none

    type(image),intent(inout)    :: img  ! Image
    type(photon),intent(in) :: p
    real(dp),intent(in) :: x_image, y_image
    integer,intent(in) :: im, inu, io
    integer :: ix,iy,ir ! Bins
    real(dp) :: transmission

    if(img%n_stokes == 4) then
       img%tmp_stokes = [p%s%i,p%s%q,p%s%u,p%s%v]
    else
       img%tmp_stokes = [p%s%i]
    end if

    if(inu >= 1 .and. inu <= img%n_nu) then
       if(img%compute_image) then
          call find_image_bin(img,x_image,y_image,ix,iy)
          if(ix >= 1 .and. ix <= img%n_x) then
             if(iy >= 1 .and. iy <= img%n_y) then
                img%img(inu,ix,iy,im,io,:) = img%img(inu,ix,iy,im,io,:) + img%tmp_stokes * p%energy * transmission
                if(img%uncertainties) then
                   img%img2(inu,ix,iy,im,io,:) = img%img2(inu,ix,iy,im,io,:) + (img%tmp_stokes * p%energy * transmission) ** 2._dp
                   img%imgn(inu,ix,iy,im,io,:) = img%imgn(inu,ix,iy,im,io,:) + 1._dp
                end if
             end if
          end if
       end if
       if(img%compute_sed) then
          call find_sed_bin(img,x_image,y_image,ir)
          if(ir >= 1 .and. ir <= img%n_ap) then
             img%sed(inu,ir,im,io,:) = img%sed(inu,ir,im,io,:) + img%tmp_stokes * p%energy * transmission
             if(img%uncertainties) then
                img%sed2(inu,ir,im,io,:) = img%sed2(inu,ir,im,io,:) + (img%tmp_stokes * p%energy * transmission) ** 2._dp
                img%sedn(inu,ir,im,io,:) = img%sedn(inu,ir,im,io,:) + 1._dp
             end if
          end if
       end if
    end if

  end subroutine image_bin_single

  subroutine image_bin_raytraced(img,p,x_image,y_image,im,spectrum)

    implicit none

    type(image),intent(inout)    :: img  ! Image
    type(photon),intent(in) :: p
    real(dp),intent(in) :: x_image, y_image
    integer,intent(in) :: im ! sub-image to bin into
    real(dp),intent(in) :: spectrum(:)
    integer :: ix,iy,ir,iw,io ! Bins
    integer :: iorig

    if(img%compute_image.and..not.allocated(img%img)) call error('image_bin_raytraced','Image not allocated')
    if(img%compute_sed.and..not.allocated(img%sed)) call error('image_bin_raytraced','SED not allocated')

    if(img%use_filters) call error("image_bin_raytraced", "filter convolution cannot be used with raytracing")

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
       iorig = orig(p)
       io = ((iorig - mod(iorig,2)) * img%n_sources + (iorig - mod(iorig+1,2) - 1) * img%n_dust) / 2
       if(mod(iorig,2)==0) then
          io = io + p%dust_id
       else
          io = io + p%source_id
       end if
    else if(trim(img%track_origin) == 'scatterings') then
       ! The first slice stores photons with no scattering, then slice 2, 3,
       ! etc. store the first scatterings, then a slice is added for the
       ! remaining flux.
       if(p%n_scat > img%track_n_scat) then
          io = img%track_n_scat + 2
       else
          io = p%n_scat + 1
       end if
       if(p%reprocessed) io = io + (img%track_n_scat + 2)
    else if(trim(img%track_origin) == 'basic') then
       io = orig(p)
    else
       io = 1
    end if

    if(img%compute_image) then
       call find_image_bin(img,x_image,y_image,ix,iy)
       if(ix >= 1 .and. ix <= img%n_x) then
          if(iy >= 1 .and. iy <= img%n_y) then
             do iw=1,img%n_nu
                img%img(iw,ix,iy,im,io,1) = img%img(iw,ix,iy,im,io,1) + spectrum(iw)
                if(img%uncertainties) then
                   img%img2(iw,ix,iy,im,io,1) = img%img2(iw,ix,iy,im,io,1) + spectrum(iw)**2._dp
                   img%imgn(iw,ix,iy,im,io,1) = img%imgn(iw,ix,iy,im,io,1) + 1._dp
                end if
             end do
          end if
       end if
    end if

    if(img%compute_sed) then
       call find_sed_bin(img,x_image,y_image,ir)
       if(ir >= 1 .and. ir <= img%n_ap) then
          do iw=1,img%n_nu
             img%sed(iw,ir,im,io,1) = img%sed(iw,ir,im,io,1) + spectrum(iw)
             if(img%uncertainties) then
                img%sed2(iw,ir,im,io,1) = img%sed2(iw,ir,im,io,1) + spectrum(iw)**2._dp
                img%sedn(iw,ir,im,io,1) = img%sedn(iw,ir,im,io,1) + 1._dp
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
    !
    ! If using filters, we want to keep the flux in F_nu * dnu for now since
    ! the filter already included the correct normalization.

    if(img%use_filters) then
      dnunorm = 1._dp
    else
      dnunorm = (img%nu_max / img%nu_min) ** (+0.5_dp / real(img%n_nu, dp)) &
           &  - (img%nu_max / img%nu_min) ** (-0.5_dp / real(img%n_nu, dp))
    end if

    if(img%compute_sed) then

       write(*,'(" [image_write] writing out SEDs")')

       allocate(cube5d(img%n_nu, img%n_ap, img%n_view,img%n_orig,img%n_stokes))
       if(img%uncertainties) allocate(cube5de(img%n_nu, img%n_ap, img%n_view,img%n_orig,img%n_stokes))

       cube5d = img%sed

       if(img%uncertainties) then
          where(img%sedn > 1)
             cube5de = sqrt((img%sed2 + (img%sed)**2 / img%sedn) / (img%sedn - 1)) * sqrt(img%sedn)
          elsewhere
             cube5de = 0._dp
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
          if(.not.img%use_filters) then
             call mp_write_keyword(group, 'seds','numin',img%nu_min)
             call mp_write_keyword(group, 'seds','numax',img%nu_max)
          end if
       end if
       call mp_write_keyword(group, 'seds','apmin',img%ap_min)
       call mp_write_keyword(group, 'seds','apmax',img%ap_max)

       call mp_write_keyword(group, 'seds', 'track_origin', img%track_origin)
       if(trim(img%track_origin) == 'detailed') then
          call mp_write_keyword(group, 'seds', 'n_sources', img%n_sources)
          call mp_write_keyword(group, 'seds', 'n_dust', img%n_dust)
       else if(trim(img%track_origin) == 'scatterings') then
          call mp_write_keyword(group, 'seds', 'track_n_scat', img%track_n_scat)
       end if

    end if

    if(img%compute_image) then

       write(*,'(" [image_write] writing out images")')

       allocate(cube6d(img%n_nu, img%n_x, img%n_y, img%n_view, img%n_orig, img%n_stokes))
       if(img%uncertainties) allocate(cube6de(img%n_nu, img%n_x, img%n_y, img%n_view, img%n_orig, img%n_stokes))

       cube6d = img%img

       if(img%uncertainties) then
          where(img%imgn > 1)
             cube6de = sqrt((img%img2 + (img%img)**2 / img%imgn) / (img%imgn - 1)) * sqrt(img%imgn)
          elsewhere
             cube6de = 0._dp
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
          if(.not.img%use_filters) then
             call mp_write_keyword(group, 'images','numin',img%nu_min)
             call mp_write_keyword(group, 'images','numax',img%nu_max)
          end if
       end if
       call mp_write_keyword(group, 'images','xmin',img%x_min)
       call mp_write_keyword(group, 'images','xmax',img%x_max)
       call mp_write_keyword(group, 'images','ymin',img%y_min)
       call mp_write_keyword(group, 'images','ymax',img%y_max)

       call mp_write_keyword(group, 'images', 'track_origin', img%track_origin)
       if(trim(img%track_origin) == 'detailed') then
          call mp_write_keyword(group, 'images', 'n_sources', img%n_sources)
          call mp_write_keyword(group, 'images', 'n_dust', img%n_dust)
       else if(trim(img%track_origin) == 'scatterings') then
          call mp_write_keyword(group, 'images', 'track_n_scat', img%track_n_scat)
       end if

    end if

    if(img%use_filters) then
      call mp_write_keyword(group, '.', 'use_filters',img%use_filters)
      call mp_write_keyword(group, '.', 'n_filt',size(img%filters))
      call mp_write_array(group, 'filt_nu0', img%filters(:)%nu0)
    end if

    if(img%use_exact_nu) then
       call mp_table_write_header(group, 'frequencies',img%n_nu,1,(/'nu'/),(/1/),(/h5t_ieee_f64le/))
       call mp_table_write_column(group, 'frequencies','nu',img%nu)
    end if

    write(*,'(" [image_write] done")')

  end subroutine image_write

end module type_image
