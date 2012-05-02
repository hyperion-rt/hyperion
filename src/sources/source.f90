module sources

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use type_source
  use type_photon
  use dust_main
  use type_vector3d
  use settings

  use grid_geometry, only : place_in_cell

  implicit none
  save

  private
  public :: setup_sources
  public :: add_source
  public :: find_nearest_source
  public :: emit
  public :: emit_peeloff

  integer,public :: n_sources = 0
  ! number of sources

  type(source),allocatable,dimension(:),public :: s
  ! the sources

  type(pdf_discrete_dp),public :: luminosity
  ! relative luminosity of all the sources

  real(dp),public :: energy_total   = 0._dp
  real(dp),public :: energy_current = 0._dp

  integer :: is
  ! loop variable

  integer(idp) :: photon_counter = 0

  logical :: any_intersect = .false.

  logical, public :: sample_sources_evenly = .false.

contains

  subroutine setup_sources(group)

    implicit none

    integer(hid_t),intent(in) :: group
    character(len=30),allocatable :: source_names(:)
    integer(hid_t) :: g_source

    if(main_process()) write(*,'(" [sources] setting up sources")')

    call mp_list_groups(group, '.', source_names)

    n_sources = size(source_names)
    allocate(s(n_sources))

    do is=1,n_sources
       g_source = mp_open_group(group, source_names(is))
       call source_read(g_source, s(is))
       call mp_close_group(g_source)
       if(s(is)%intersect) any_intersect = .true.
    end do

    ! flux = flux / nu
    ! energy = h*nu
    ! flux should represent number of photons or at least the relative number. Normalise such that integral = number of photons
    ! then flux * h*nu gives total energy without normalisation, luminosity over that gives normalization factor.

    if(n_sources > 0) then
       call set_pdf(luminosity,s(:)%luminosity)
       energy_total = sum(s(:)%luminosity)
    else
       energy_total = 0._dp
    end if

  end subroutine setup_sources

  subroutine add_source(new_source)
    implicit none
    type(source),intent(in) :: new_source
    type(source),allocatable :: tmp_s(:)
    allocate(tmp_s(n_sources))
    tmp_s = s
    deallocate(s)
    allocate(s(n_sources + 1))
    s(1:n_sources) = tmp_s(1:n_sources)
    s(n_sources+1) = new_source
    n_sources = n_sources + 1
    if(new_source%intersect) any_intersect = .true.
    call set_pdf(luminosity,s(:)%luminosity)
    energy_total = sum(s(:)%luminosity)
  end subroutine add_source

  subroutine emit(p,reemit,reemit_id,reemit_energy,inu)

    implicit none

    type(photon),intent(out) :: p
    ! the emitted photon

    logical,intent(in),optional :: reemit
    integer,intent(in),optional :: reemit_id
    real(dp),intent(in),optional :: reemit_energy
    integer,intent(in),optional :: inu
    logical :: do_reemit
    real(dp) :: xi

    if(n_sources==0) call error("emit", "no sources to emit from")

    ! --- First, decide which source to sample from --- !

    if(n_sources > 1) then
       if(sample_sources_evenly) then
          call random(xi)
          p%source_id = int(xi*n_sources) + 1
       else
          p%source_id = sample_pdf(luminosity)
       end if
    else
       p%source_id = 1
    end if

    if(present(reemit)) then
       do_reemit = reemit
    else
       do_reemit = .false.
    end if

    if(do_reemit) then
       if(present(reemit_id)) then
          p%source_id = reemit_id
       else
          call error("emit", "reemit_id is missing")
       end if
    end if

    ! --- Then emit the photon --- !

    if(present(inu)) then
       call source_emit(s(p%source_id),p,nu=frequencies(inu))
       p%inu = inu
    else
       call source_emit(s(p%source_id),p)
    end if

    call angle3d_to_vector3d(p%a,p%v)

    if(do_reemit) then
       if(present(reemit_energy)) then
          p%energy = reemit_energy
       else
          call error("emit", "reemit_energy is missing")
       end if
    else
       if(present(inu)) p%energy = p%energy * energy_total
       if(sample_sources_evenly) p%energy = p%energy * luminosity%pdf(p%source_id) * n_sources
       energy_current = energy_current + p%energy
    end if

    call prepare_photon(p)
    call update_optconsts(p)

    p%emiss_type = s(p%source_id)%freq_type

    photon_counter = photon_counter + 1
    p%id = photon_counter

    p%last = 'sr'

    call place_in_cell(p)
    if(p%killed) call error("emit", "photon was not emitted inside a cell - this usually indicates that a source is not inside the grid")

  end subroutine emit

  subroutine emit_peeloff(p,a_req)
    implicit none
    type(photon),intent(inout) :: p ! the photon to peeloff
    type(angle3d_dp),intent(in) :: a_req ! requested angle
    call source_emit_peeloff(s(p%source_id),p,a_req)
    call angle3d_to_vector3d(p%a,p%v)
  end subroutine emit_peeloff

  logical function intersects_sources(r1,r2)

    type(vector3d_dp),intent(in) :: r1,r2
    integer :: is

    do is=1,n_sources
       if(source_intersect(s(is),r1,r2)) then
          intersects_sources = .true.
          return
       end if
    end do

    intersects_sources = .false.
    return

  end function intersects_sources

  subroutine find_nearest_source(r,v,nearest_distance,nearest_id)

    type(vector3d_dp),intent(in) :: r,v
    integer :: is
    real(dp),intent(out) :: nearest_distance
    integer,intent(out) :: nearest_id

    nearest_id = 0
    nearest_distance = infinity_dp()

    if(.not.any_intersect) return

    do is=1,n_sources
       if(s(is)%intersect) then
          if(source_distance(s(is),r,v) < nearest_distance) then
             nearest_distance = source_distance(s(is),r,v)
             nearest_id = is
          end if
       end if
    end do

  end subroutine find_nearest_source

end module sources

