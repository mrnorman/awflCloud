
KOKKOS_PATH = /home/imn/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = driver.cpp

default: main

CXX := mpic++
# CFLAGS := -DARRAY_DEBUG -O1 -g -I${PNETCDF_PATH}/include
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

KOKKOS_DEVICES = "Serial"

include ${KOKKOS_PATH}/Makefile.kokkos

main: $(KOKKOS_LINK_DEPENDS) $(KOKKOS_CPP_DEPENDS) driver.cpp
	$(CXX) ${CFLAGS} $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) driver.cpp -o cloudFV $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat cloudFV

realclean: clean
	rm -f KokkosCore_config.h KokkosCore_config.tmp libkokkos.a
