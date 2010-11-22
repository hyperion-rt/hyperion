# Info

AC_INIT(rt,version-1.0)

# Docstring for options

AC_ARG_ENABLE(profiling, [  --enable-profiling      Enables the profiling flags])
AC_ARG_ENABLE(debug, [  --enable-debug          Enables the full debugging flags])

# Profiling flags

if test "$enable_profiling" == "yes"
then
  profiling_ifort='-g -save-temps'
  profiling_f95=''
  profiling_gfortran=''
else
  profiling_ifort=''
  profiling_f95=''
  profiling_gfortran=''
fi

# Debugging flags

if test "$enable_debug" == "yes"
then
  debug_f95='-C=all -gline -nan'
  debug_ifort='-check all -warn all -warn nodec -warn interfaces -gen_interfaces -traceback -g'
  debug_gfortran='-Wall -fbounds-check -g'
else
  debug_f95=''
  debug_ifort=''
  debug_gfortran=''
fi

extra_f95="-Isrc/modules -mdir src/modules -f2003 -colour -ieee=full $debug_f95 $profiling_f95"
extra_ifort="-module src/modules -heap-arrays -L/usr/lib -lgcc -error-limit 1 $debug_ifort $profiling_ifort"
extra_gfortran="-Jsrc/modules -Isrc/modules -g -ffree-line-length-none $debug_gfortran $profiling_gfortran"

# SERIAL CODE

AC_CHECK_PROG(exists_hdf, [h5fc], yes, no)

if test "$exists_hdf" == no
then
  AC_ERROR("The HDF5 library with support for Fortran needs to be installed")
else
  AC_SUBST(serial_fc, 'h5fc')
fi

# Find out what compiler is buried under the HDF5 wrapper
hdf_compiler=`h5fc -show | awk {'print $1'}`

if test "$hdf_compiler" == f95
then
  AC_SUBST(posix_module_serial, posix_nag)
  AC_SUBST(extra_serial, "$extra_f95")
elif test "$hdf_compiler" == ifort
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_ifort")
elif test "$hdf_compiler" == gfortran
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_gfortran")
fi

# MPI CODE

# First find out if we can use h5pfc
AC_CHECK_PROGS(parallel_compiler, [mpif90 h5pfc])

if test "$parallel_compiler" == h5pfc
then
  mpi_compiler=`h5pfc -show | awk {'print $1'}`
  if test "$mpi_compiler" == mpif90
  then
    mpi_compiler=`mpif90 -show | awk {'print $1'}`
  fi
  AC_SUBST(mpi_fc, 'h5pfc')
elif test "$parallel_compiler" == mpif90
then
  mpi_compiler=`mpif90 -show | awk {'print $1'}`
  AC_SUBST(mpi_fc, 'HDF5_FLINKER=mpif90 HDF5_FC=mpif90 h5fc')
else
  AC_ERROR(["h5pfc or mpif90 need to be present"])
fi

if test "$mpi_compiler" == f95
then
  AC_SUBST(posix_module_mpi, posix_nag)
  AC_SUBST(extra_mpi, "$extra_f95")
elif test "$mpi_compiler" == ifort
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_ifort")
elif test "$mpi_compiler" == gfortran
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_gfortran")
fi

# Output final Makefile

AC_OUTPUT(Makefile)