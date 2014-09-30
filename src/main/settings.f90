module settings

  use core_lib

  implicit none
  save

  integer(idp) :: n_inter_max = 0
  logical :: n_inter_max_warn = .true.
  integer(idp) :: n_mrw_max = 0
  logical :: n_mrw_max_warn = .true.
  integer(idp) :: n_reabs_max = 0
  logical :: n_reabs_max_warn = .true.

  integer(idp),public :: n_stats
  integer,public :: n_initial_iter
  integer(idp),public :: n_initial_photons
  integer(idp),public :: n_last_photons = 0
  integer(idp),public :: n_last_photons_sources = 0
  integer(idp),public :: n_last_photons_dust = 0
  integer(idp),public :: n_raytracing_photons_sources = 0
  integer(idp),public :: n_raytracing_photons_dust = 0
  logical,public :: use_raytracing, use_mrw, use_pda
  logical, public :: kill_on_absorb, kill_on_scatter
  real(dp),public :: mrw_gamma
  logical, public :: forced_first_scattering

  logical :: use_exact_nu = .false.
  real(dp),allocatable :: frequencies(:)

  integer :: physics_io_type

  character(len=4) :: output_density
  character(len=4) :: output_density_diff
  character(len=4) :: output_specific_energy
  character(len=4) :: output_n_photons

  logical :: check_convergence = .false.
  real(dp) :: convergence_absolute = 0._dp
  real(dp) :: convergence_relative = 0._dp
  real(dp) :: convergence_percentile = 100._dp

  logical :: enforce_energy_range

  real(dp) :: propagation_check_frequency

  character(len=20) :: specific_energy_type

end module settings
