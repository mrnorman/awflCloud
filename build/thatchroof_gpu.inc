
CXX = nvcc -ccbin /usr/bin/mpic++ --expt-relaxed-constexpr -expt-extended-lambda -x cu --std=c++14 --use_fast_math -expt-extended-lambda -arch=sm_35 -D__USE_CUDA__ -I/usr/local/cuda-10.1/include
LINK = nvcc -ccbin /usr/bin/mpic++
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

