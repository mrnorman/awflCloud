
CC := mpic++
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} driver.cpp -o cloudFV ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat cloudFV

