include mach.inc

main: $(EXE)

$(EXE): $(OBJ) YAKL.o
	$(LINK) $(OBJ) $(LIB) -o $(EXE) $(LDFLAGS) YAKL.o

%.o:%.cpp *.h
	$(CXX) $(CXXFLAGS) -I./cub -c $< -o $(notdir $@)

clean:
	rm -f *.gch *.o *.dat cloudFV

realclean: clean

distclean: clean


