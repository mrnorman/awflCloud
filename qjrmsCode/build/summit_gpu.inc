
PNETCDF_PATH = ${OLCF_PARALLEL_NETCDF_PATH}

CXX = nvcc -ccbin mpic++ -expt-extended-lambda -x cu --std=c++14 --use_fast_math -arch=sm_70 -D__USE_CUDA__
LINK = nvcc -ccbin mpic++
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include -I/ccs/home/imn/cub
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

