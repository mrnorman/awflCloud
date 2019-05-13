
CC := mpic++
# CC := /usr/local/cuda-10.1/bin/nvcc --expt-extended-lambda -x cu -ccbin /usr/bin/mpic++
# CFLAGS := -DARRAY_DEBUG -O0 -g -I${PNETCDF_PATH}/include
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} driver.cpp -o cloudFV ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat cloudFV
