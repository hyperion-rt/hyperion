module mpi_io

  use core_lib, only : mp_read_keyword => hdf5_read_keyword, &
       &               mp_table_read_column_auto => hdf5_table_read_column_auto, &
       &               mp_read_array_auto => hdf5_read_array_auto

  implicit none
  save

  private
  public :: mp_read_keyword
  public :: mp_table_read_column_auto
  public :: mp_read_array_auto

end module mpi_io

