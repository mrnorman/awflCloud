include mach.inc

main: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(LDFLAGS)

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) *.h
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)

clean:
	rm -f *.gch *.o *.dat cloudFV

realclean: clean
	rm -f KokkosCore_config.h KokkosCore_config.tmp libkokkos.a

distclean: clean
	rm -f KokkosCore_config.h KokkosCore_config.tmp libkokkos.a


