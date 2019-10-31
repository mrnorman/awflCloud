include mach.inc

main: $(EXE)

$(EXE): $(OBJ)
	$(LINK) $(OBJ) $(LIB) -o $(EXE) $(LDFLAGS)

%.o:%.cpp *.h
	$(CXX) $(CXXFLAGS) -I./cub -c $< -o $(notdir $@)

clean:
	rm -f *.gch *.o *.dat cloudFV

realclean: clean

distclean: clean


