module counters

  use core_lib, only : idp

  implicit none
  save

  integer(idp) :: killed_photons_geo = 0
  integer(idp) :: killed_photons_int = 0
  integer(idp) :: killed_photons_reabs = 0

end module counters

