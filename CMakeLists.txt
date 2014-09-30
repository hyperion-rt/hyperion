# CMake project file for Hyperion.
cmake_minimum_required(VERSION 2.8)
project(hyperion)
enable_language(Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Make sure that the default is a Release.
if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif()

# Options setup.
option(ENABLE_PROFILING "Enable profiling." OFF)
option(ENABLE_STATIC "Enable static compilation." OFF)

# Find the HDF5 packages.
find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL HL)
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
set(MANDATORY_LIBRARIES ${HDF5_LIBRARIES})
message(STATUS "HDF5 is parallel: ${HDF5_IS_PARALLEL}")

# MPI setup.
find_package(MPI REQUIRED)
if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Could not detect the MPI Fortran environment, aborting.")
endif()
message(STATUS "MPI include dirs: ${MPI_Fortran_INCLUDE_PATH}")
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Compiler detection.
message(STATUS "Compiler id: ${CMAKE_Fortran_COMPILER_ID}")

# Flags setup.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # There are no special flags for profiling or static compilation when using gfortran.
    # Base flags used in Debug and Release.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -ffree-line-length-none")
    # Debug-specific flags.
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -fbounds-check")
    # POSIX module setup.
    set(HYPERION_posix_module_serial "posix_default.f90")
    set(HYPERION_posix_module_mpi "posix_default.f90")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # ifort.
    set(HYPERION_Fortran_FLAGS "-heap-arrays -error-limit 1 -fp-model strict -diag-disable 8291")
    if(ENABLE_PROFILING)
        # Intel profiling flags.
        set(HYPERION_Fortran_FLAGS "${HYPERION_Fortran_FLAGS} -g -save-temps")
    endif()
    if(ENABLE_STATIC)
        # Intel static compilation flags.
        set(HYPERION_Fortran_FLAGS "${HYPERION_Fortran_FLAGS} --static")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HYPERION_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -warn all -warn nodec -warn interfaces -gen_interfaces -traceback -g")
    # POSIX module setup.
    set(HYPERION_posix_module_serial "posix_default.f90")
    set(HYPERION_posix_module_mpi "posix_default.f90")
else()
    message(FATAL_ERROR "Could not determine the compiler, or compiler is unsupported.")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Release compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

set(COMMON_FILES
    fortranlib/src/base_types.f90
    fortranlib/src/lib_messages.f90
    fortranlib/src/lib_io.f90
    fortranlib/src/lib_version.f90
    fortranlib/src/lib_algebra.f90
    fortranlib/src/lib_conf.f90
    fortranlib/src/lib_array.f90
    fortranlib/src/lib_statistics.f90
    fortranlib/src/lib_random.f90
    fortranlib/src/lib_hdf5.f90
    fortranlib/src/lib_constants.f90
    fortranlib/src/type_angle3d.f90
    fortranlib/src/type_vector3d.f90
    fortranlib/src/type_stokes.f90
    fortranlib/src/type_pdf.f90
)

set(OBJECTS1_FILES
    src/core/type_photon.f90
    src/main/settings.f90
    src/main/performance.f90
    src/images/image_type.f90
    src/dust/dust_type_4elem.f90
    src/dust/dust.f90
)

set(OBJECTS2_FILES
    src/grid/grid_physics_3d.f90
    src/grid/grid_monochromatic.f90
    src/dust/dust_interact.f90
    src/sources/source_type.f90
    src/sources/source.f90
    src/grid/grid_propagate_3d.f90
    src/grid/grid_mrw_3d.f90
    src/grid/grid_pda_3d.f90
    src/grid/grid_generic.f90
    src/images/images_binned.f90
    src/images/images_peeled.f90
    src/main/setup_rt.f90
)

set(OBJECTS3_FILES
    src/main/iter_lucy.f90
    src/main/iter_final.f90
    src/main/iter_final_mono.f90
    src/main/iter_raytracing.f90
    src/main/main.f90
)

set(OBJECTS_CAR_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_cartesian_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_cartesian_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_cartesian_3d.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_CAR_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_cartesian_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_cartesian_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_cartesian_3d.f90
    src/mpi/mpi_routines.f90
)

set(OBJECTS_AMR_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_amr.f90
    src/core/type_cell_id_amr.f90
    src/grid/grid_io_amr.f90
    src/grid/grid_geometry_amr.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_AMR_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_amr.f90
    src/core/type_cell_id_amr.f90
    src/grid/grid_io_amr.f90
    src/grid/grid_geometry_amr.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines.f90
)

set(OBJECTS_OCT_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_octree.f90
    src/core/type_cell_id_octree.f90
    src/grid/grid_io_1d.f90
    src/grid/grid_geometry_octree.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_OCT_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_octree.f90
    src/core/type_cell_id_octree.f90
    src/grid/grid_io_1d.f90
    src/grid/grid_geometry_octree.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines.f90
)

set(OBJECTS_VOR_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/extern/kdtree2/kdtree2.f90
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_voronoi.f90
    src/core/type_cell_id_octree.f90
    src/grid/grid_io_1d.f90
    src/grid/grid_geometry_voronoi.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_VOR_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/extern/kdtree2/kdtree2.f90
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_voronoi.f90
    src/core/type_cell_id_octree.f90
    src/grid/grid_io_1d.f90
    src/grid/grid_geometry_voronoi.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_disabled.f90
    src/mpi/mpi_routines.f90
)

set(OBJECTS_SPH_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_spherical_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_spherical_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_spherical_3d.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_SPH_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_spherical_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_spherical_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_spherical_3d.f90
    src/mpi/mpi_routines.f90
)

set(OBJECTS_CYL_SER_FILES
    src/mpi/mpi_core_nompi.f90
    fortranlib/src/${HYPERION_posix_module_serial}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io_nompi.f90
    src/grid/type_grid_cylindrical_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_cylindrical_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_cylindrical_3d.f90
    src/mpi/mpi_routines_nompi.f90
)

set(OBJECTS_CYL_MPI_FILES
    src/mpi/mpi_core.f90
    fortranlib/src/${HYPERION_posix_module_mpi}
    src/core/core_lib.f90
    src/main/counters.f90
    src/mpi/mpi_io.f90
    src/grid/type_grid_cylindrical_3d.f90
    src/core/type_cell_id_3d.f90
    src/grid/grid_io.f90
    src/grid/grid_geometry_cylindrical_3d.f90
    src/grid/grid_geometry_common_3d.f90
    src/grid/grid_pda_cylindrical_3d.f90
    src/mpi/mpi_routines.f90
)

macro(HYPERION_ADD_SERIAL_EXECUTABLE name files)
    set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/hyperion_${name}_mods/") 
    add_executable("hyperion_${name}" ${COMMON_FILES} ${OBJECTS1_FILES} ${OBJECTS2_FILES} ${OBJECTS3_FILES} ${files})
    target_link_libraries("hyperion_${name}" ${MANDATORY_LIBRARIES})
    install(TARGETS "hyperion_${name}" RUNTIME DESTINATION bin)
endmacro()

macro(HYPERION_ADD_MPI_EXECUTABLE name files)
    set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/hyperion_${name}_mpi_mods/") 
    add_executable("hyperion_${name}_mpi" ${COMMON_FILES} ${OBJECTS1_FILES} ${OBJECTS2_FILES} ${OBJECTS3_FILES} ${files})
    target_link_libraries("hyperion_${name}_mpi" ${MANDATORY_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    install(TARGETS "hyperion_${name}_mpi" RUNTIME DESTINATION bin)
endmacro()

HYPERION_ADD_SERIAL_EXECUTABLE(car "${OBJECTS_CAR_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(car "${OBJECTS_CAR_MPI_FILES}")
HYPERION_ADD_SERIAL_EXECUTABLE(amr "${OBJECTS_AMR_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(amr "${OBJECTS_AMR_MPI_FILES}")
HYPERION_ADD_SERIAL_EXECUTABLE(oct "${OBJECTS_OCT_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(oct "${OBJECTS_OCT_MPI_FILES}")
HYPERION_ADD_SERIAL_EXECUTABLE(vor "${OBJECTS_VOR_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(vor "${OBJECTS_VOR_MPI_FILES}")
HYPERION_ADD_SERIAL_EXECUTABLE(sph "${OBJECTS_SPH_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(sph "${OBJECTS_SPH_MPI_FILES}")
HYPERION_ADD_SERIAL_EXECUTABLE(cyl "${OBJECTS_CYL_SER_FILES}")
HYPERION_ADD_MPI_EXECUTABLE(cyl "${OBJECTS_CYL_MPI_FILES}")
