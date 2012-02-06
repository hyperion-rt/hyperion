# Info

AC_INIT(rt,version-1.0)

# Docstring for options

AC_ARG_ENABLE(profiling, [  --enable-profiling      Enables the profiling flags])
AC_ARG_ENABLE(debug, [  --enable-debug          Enables the full debugging flags])
AC_ARG_ENABLE(static, [  --enable-static         Enables static compilation])

# Static flags

if test "$enable_static" == "yes"
then
  static_ifort='--static'
  static_nag=''
  static_gfortran=''
  static_g95=''
  static_pgi='-Bstatic'
else
  static_ifort=''
  static_nag=''
  static_gfortran=''
  static_g95=''
  static_pgi=''
fi

# Profiling flags

if test "$enable_profiling" == "yes"
then
  profiling_ifort='-g -save-temps'
  profiling_nag=''
  profiling_gfortran=''
  profiling_g95=''
  profiling_pgi=''
else
  profiling_ifort=''
  profiling_nag=''
  profiling_gfortran=''
  profiling_g95=''
  profiling_pgi=''
fi

# Debugging flags

if test "$enable_debug" == "yes"
then
  debug_nag='-C=all -gline -nan'
  debug_ifort='-check all -warn all -warn nodec -warn interfaces -gen_interfaces -traceback -g'
  debug_gfortran='-Wall -fbounds-check -g'
  debug_g95='-Wall -fbounds-check -g'
  debug_pgi='-Mbounds -traceback -g'
else
  debug_nag=''
  debug_ifort=''
  debug_gfortran=''
  debug_g95=''
  debug_pgi=''
fi

extra_nag="-Isrc/modules -mdir src/modules -f2003 -colour -ieee=full $debug_nag $profiling_nag $static_nag"
extra_ifort="-module src/modules -heap-arrays -L/usr/lib -lgcc -error-limit 1 $debug_ifort $profiling_ifort $static_ifort"
extra_gfortran="-Jsrc/modules -Isrc/modules -g -ffree-line-length-none $debug_gfortran $profiling_gfortran $static_gfortran"
extra_g95="-fmod=src/modules -Isrc/modules -g -ffree-line-length-huge $debug_g95 $profiling_g95 $static_g95"
extra_pgi="$debug_pgi $profiling_pgi $static_pgi"

# SERIAL CODE

AC_CHECK_PROGS(hdf_compiler, [h5fc h5pfc], none)

if test "$hdf_compiler" == none
then
  AC_ERROR("The HDF5 library with support for Fortran needs to be installed")
else
  AC_SUBST(serial_fc, "$hdf_compiler")
fi

# Find out what compiler is buried under the HDF5 wrapper
if test "$hdf_compiler" == "h5fc"
then
  hdf_compiler=`h5fc -show | awk {'print $1'}`
else
  hdf_compiler=`h5pfc -show | awk {'print $1'}`
  if test "$hdf_compiler" == "mpif90"
  then
    hdf_compiler=`mpif90 -show | awk {'print $1'}`
  fi
fi

# Remove any leading path from the compiler
hdf_compiler=`basename $hdf_compiler`

# gfortran sometimes includes version number, so check whether this is the case
hdf_compiler_short=`echo $hdf_compiler | cut -c1-8`
if test "$hdf_compiler_short" == gfortran
then
  hdf_compiler=gfortran
fi

if test "$hdf_compiler" == f95
then
  AC_SUBST(posix_module_serial, posix_nag)
  AC_SUBST(extra_serial, "$extra_nag")
elif test "$hdf_compiler" == nagfor
then
  AC_SUBST(posix_module_serial, posix_nag)
  AC_SUBST(extra_serial, "$extra_nag")
elif test "$hdf_compiler" == pgf90
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$hdf_compiler" == pgf95
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$hdf_compiler" == pgfortran
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_pgi")
elif test "$hdf_compiler" == ifort
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_ifort")
elif test "$hdf_compiler" == gfortran
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_gfortran")
elif test "$hdf_compiler" == g95
then
  AC_SUBST(posix_module_serial, posix_default)
  AC_SUBST(extra_serial, "$extra_g95")
fi

# MPI CODE

# First find out if we can use h5pfc
AC_CHECK_PROGS(parallel_compiler, [h5pfc mpif90])

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

# Remove any leading path from the compiler
mpi_compiler=`basename $mpi_compiler`

# gfortran sometimes includes version number, so check whether this is the case
mpi_compiler_short=`echo $mpi_compiler | cut -c1-8`
if test "$mpi_compiler_short" == gfortran
then
  mpi_compiler=gfortran
fi

if test "$mpi_compiler" == f95
then
  AC_SUBST(posix_module_mpi, posix_nag)
  AC_SUBST(extra_mpi, "$extra_nag")
elif test "$mpi_compiler" == nagfor
then
  AC_SUBST(posix_module_mpi, posix_nag)
  AC_SUBST(extra_mpi, "$extra_nag")
elif test "$mpi_compiler" == pgf90
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_pgi")
elif test "$mpi_compiler" == pgf95
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_pgi")
elif test "$mpi_compiler" == pgfortran
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_pgi")
elif test "$mpi_compiler" == ifort
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_ifort")
elif test "$mpi_compiler" == gfortran
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_gfortran")
elif test "$mpi_compiler" == g95
then
  AC_SUBST(posix_module_mpi, posix_default)
  AC_SUBST(extra_mpi, "$extra_g95")
fi

# Output final Makefile

AC_OUTPUT(Makefile)