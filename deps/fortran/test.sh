#!/bin/bash

h5fc test_hdf5.f90 -o test_hdf5
./test_hdf5

mpif90 test_mpi.f90 -o test_mpi
./test_mpi

HDF5_FLINKER=mpif90 HDF5_FC=mpif90 h5fc test_hdf5_mpi.f90 -o test_hdf5_mpi
./test_hdf5_mpi
