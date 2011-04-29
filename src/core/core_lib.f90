module core_lib

  use posix, only : microsleep

  ! Configuration file I/O
  use lib_conf
  use lib_io
  use lib_version

  ! Error/message handling
  use lib_messages

  ! Basic numeric types
  use base_types

  ! Maths
  use lib_constants
  use lib_algebra
  use lib_statistics
  use type_pdf
  use type_angle3d
  use type_vector3d
  use type_stokes
  use lib_array
  use lib_random

  ! External
  use lib_hdf5

end module core_lib
