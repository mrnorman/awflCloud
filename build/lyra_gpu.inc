
SRC = driver.cpp
EXE = cloudFV

PNETCDF_PATH=/ccs/home/imn/parallel-netcdf-1.11.2_clang

default: main

CXX = mpic++ -std=c++14
LINK = mpic++
CXXFLAGS = -O3 -ffast-math -I${PNETCDF_PATH}/include -D__USE_HIP__
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

