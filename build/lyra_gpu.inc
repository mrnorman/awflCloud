
SRC = driver.cpp
EXE = cloudFV

PNETCDF_PATH=/ccs/home/imn/parallel-netcdf-1.11.2_clang

default: main

CXX = hipcc -std=c++14
LINK = hipcc
CXXFLAGS = --amdgpu-target=gfx900 -g -I${PNETCDF_PATH}/include -I/opt/rocm/include/hip -D__USE_HIP__ -I/opt/openmpi/4.0.0-ucx/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf -L/opt/openmpi/4.0.0-ucx/lib -lmpi

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

