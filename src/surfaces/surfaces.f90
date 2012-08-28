module surface_collection

  use core_lib
  use mpi_core
  use mpi_hdf5_io

  use type_surface, only : surface, surface_intersect, surface_distance, surface_read

  implicit none
  save

  private
  public :: setup_surfaces
  public :: register_surface
  public :: find_nearest_surface

  integer, public :: n_surfaces
  ! number of surfaces

  type(surface), allocatable, target, public :: surfaces(:)
  ! the surfaces

  integer :: is
  ! loop variable

contains

  subroutine setup_surfaces(group)

    ! Set up the surfaces described in an HDF5 group
    !
    ! Parameters
    ! ----------
    ! group : HDF5 group
    !    The group to read the surfaces from

    integer(hid_t), intent(in) :: group

    character(len=30), allocatable :: surface_names(:)
    integer(hid_t) :: g_surface

    if(main_process()) write(*, '(" [surfaces] setting up surfaces")')

    call mp_list_groups(group, '.', surface_names)

    n_surfaces = size(surface_names)
    allocate(surfaces(n_surfaces))

    do is=1, n_surfaces
       g_surface = mp_open_group(group, surface_names(is))
       call surface_read(g_surface, surfaces(is))
       call mp_close_group(g_surface)
    end do

  end subroutine setup_surfaces

  logical function intersects_surfaces(r1, r2)

    ! Finds out whether a segment intersects any surface
    !
    ! Parameters
    ! ----------
    ! r1, r2 : vector3d_dp
    !     The positions defining the extremities of the segment

    type(vector3d_dp), intent(in) :: r1, r2

    do is=1, n_surfaces
       if(surface_intersect(surfaces(is), r1, r2)) then
          intersects_surfaces = .true.
          return
       end if
    end do

    intersects_surfaces = .false.

  end function intersects_surfaces

  subroutine find_nearest_surface(r, v, nearest_distance, nearest_id)

    ! Find the first surface that is intersected
    !
    ! Parameters
    ! ----------
    ! r : vector3d_dp
    !     The current position of the photon
    ! v : vector3d_dp
    !     The direction vector of the photon
    !
    ! Returns
    ! -------
    ! nearest_distance : real(dp)
    !     The distance to the nearest surface
    ! nearest_id : integer
    !     The index of the nearest surface

    type(vector3d_dp), intent(in) :: r, v

    real(dp), intent(out) :: nearest_distance
    integer, intent(out) :: nearest_id

    nearest_id = 0
    nearest_distance = infinity_dp()

    do is=1, n_surfaces
       if(surface_distance(surfaces(is), r, v) < nearest_distance) then
          nearest_distance = surface_distance(surfaces(is), r, v)
          nearest_id = is
       end if
    end do

  end subroutine find_nearest_surface

  subroutine register_surface(srf)

    ! Register a surface with the collection
    !
    ! Parameters
    ! ----------
    ! srf : surface object
    !     The surface object to register with the surface collection

    type(surface), intent(in) :: srf

    type(surface), allocatable :: surfaces_tmp(:)

    ! Backup current surface list
    allocate(surfaces_tmp(n_surfaces))
    surfaces_tmp = surfaces

    ! Increase size of surface list
    deallocate(surfaces)
    allocate(surfaces(n_surfaces + 1))
    surfaces(1:n_surfaces) = surfaces_tmp

    ! Add new surface to list
    surfaces(n_surfaces + 1) = srf

    ! Increase surface counter
    n_surfaces = n_surfaces + 1

  end subroutine register_surface


end module surface_collection
