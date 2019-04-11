
CC := mpic++
# CFLAGS := -DARRAY_DEBUG -O1 -g -I${PNETCDF_PATH}/include
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} driver.cpp -o cloudFV ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat cloudFV
