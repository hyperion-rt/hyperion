! Verified and cleaned 13/09/08

module type_photon

  use core_lib
  use type_grid_cell

  implicit none
  save

  private
  public :: photon

  type photon

     type(vector3d_dp) :: r       ! Position
     type(vector3d_dp) :: v       ! Direction (vector)
     type(angle3d_dp)  :: a       ! Direction (angle)
     type(stokes_dp)   :: s       ! Stokes parameters
     real(dp)          :: nu      ! Frequency
     real(dp)          :: energy  ! Energy
     logical(1) :: in_cell = .false. ! Whether the photon is in a cell
     logical(1) :: on_wall = .false. ! Whether the photon is on a wall
     type(wall_id)     :: on_wall_id = no_wall

     integer(idp) :: id = 0 ! Photon ID (unique)

     logical :: killed = .false.
     logical :: reabsorbed = .false.
     integer :: reabsorbed_id = 0

     real(dp),allocatable  :: current_chi(:)
     real(dp),allocatable  :: current_albedo(:)
     real(dp),allocatable  :: current_kappa(:)

     type(grid_cell) :: icell

     character(len=2) :: last
     type(angle3d_dp) :: a_prev
     type(stokes_dp) :: s_prev
     type(vector3d_dp) :: v_prev
     logical :: last_isotropic

     ! RAYTRACING

     ! The emissivity type:
     ! 1 - spectrum from a source (use source_id)
     ! 2 - blackbody (use temperature)
     ! 3 - dust emission (use specific_energy and dust_id)
     integer :: emiss_type = 0

     integer :: source_id = 0
     type(angle3d_dp) :: source_a

     integer :: dust_id = 0
     real(dp) :: emiss_var_frac = 0.
     integer :: emiss_var_id = 0

     ! The following flag keeps track of whether the photon has been
     ! reprocessed (absorbed and re-emitted) by dust yet
     logical :: reprocessed=.false.

     ! The following flag keeps track of whether the photon has been scattered
     ! since being last emitted by dust or sources, and the total number of 
     ! scatterings
     logical :: scattered = .false.
     integer :: n_scat = 0

     integer :: inu

     integer :: face_id = 0

  end type photon

end module type_photon
