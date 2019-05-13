
# CC := mpic++
CC := nvcc --expt-extended-lambda -x cu -ccbin ${OLCF_SPECTRUM_MPI_ROOT}/bin/mpic++
# CFLAGS := -DARRAY_DEBUG -O0 -g -I${PNETCDF_PATH}/include
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} driver.cpp -o cloudFV ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat cloudFV
