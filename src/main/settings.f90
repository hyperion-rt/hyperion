module settings

  use core_lib

  implicit none
  save

  real(dp) :: minimum_temperature = 0.
  integer :: dust_sublimation_mode = 0 ! 0=no, 1=fast, 2=slow
  real(dp) :: dust_sublimation_temperature
  integer(idp) :: n_inter_max = 0
  integer(idp) :: n_mrw_max = 0
  integer(idp) :: n_reabs_max = 0

  integer(idp),public :: n_stats
  ! integer(idp),public :: n_lucy_check
  integer,public :: n_lucy_iter
  integer(idp),public :: n_lucy_photons
  integer(idp),public :: n_last_photons = 0
  integer(idp),public :: n_last_photons_sources = 0
  integer(idp),public :: n_last_photons_dust = 0
  integer(idp),public :: n_raytracing_photons_sources = 0
  integer(idp),public :: n_raytracing_photons_dust = 0
  logical,public :: use_raytracing, use_mrw, use_pda
  logical, public :: kill_on_absorb
  real(dp),public :: mrw_gamma
  logical, public :: forced_first_scattering

  logical :: use_exact_nu = .false.
  real(dp),allocatable :: frequencies(:)

  integer :: physics_io_type

  character(len=4) :: output_temperature
  character(len=4) :: output_density
  character(len=4) :: output_specific_energy_abs
  character(len=4) :: output_n_photons

  logical :: check_convergence = .false.
  real(dp) :: convergence_absolute = 0._dp
  real(dp) :: convergence_relative = 0._dp
  real(dp) :: convergence_percentile = 100._dp

end module settings
