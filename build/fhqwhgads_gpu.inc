
KOKKOS_PATH = ${HOME}/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = driver.cpp
EXE = cloudFV

default: main

CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper -x cu -ccbin /usr/bin/mpic++
LINK = ${CXX}
# CFLAGS = -DARRAY_DEBUG -O1 -g -I${PNETCDF_PATH}/include
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

KOKKOS_DEVICES = "Cuda"
KOKKOS_ARCH = "Maxwell50"
KOKKOS_CUDA_OPTIONS += "enable_lambda"

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include ${KOKKOS_PATH}/Makefile.kokkos

