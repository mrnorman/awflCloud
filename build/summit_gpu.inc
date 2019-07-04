
KOKKOS_PATH = /ccs/home/imn/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = driver.cpp
EXE = cloudFV

default: main

CXX = ./nvcc_wrapper --use_fast_math -x cu -ccbin ${OLCF_SPECTRUM_MPI_ROOT}/bin/mpic++
LINK = ${CXX}
# CFLAGS = -DARRAY_DEBUG -O1 -g -I${PNETCDF_PATH}/include
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

KOKKOS_DEVICES = "Cuda"
KOKKOS_ARCH = "Volta70"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
KOKKOS_OPTIONS = "disable_profiling"

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include ${KOKKOS_PATH}/Makefile.kokkos

