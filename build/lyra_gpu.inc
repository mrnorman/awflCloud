
SRC = driver.cpp
EXE = cloudFV

PNETCDF_PATH=/ccs/home/imn/parallel-netcdf-1.11.2_clang

default: main

CXX = hipcc -std=c++14
LINK = hipcc
CXXFLAGS = -amdgpu-target=mi25 -g -I${PNETCDF_PATH}/include -I/opt/rocm/include/hip -D__USE_HIP__ -I/ccs/home/imn/openmpi_clang/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf -L/ccs/home/imn/openmpi_clang/lib -lmpi

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

