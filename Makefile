include mach.inc

EXE = cloudFV
SOURCES = driver.cpp FileIO.cpp params.cpp YAKL.cpp Exchange.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(info $(OBJECTS))

default: main

main: $(EXE)

$(EXE): $(OBJECTS)
	$(LINK) $(OBJECTS) $(LIB) -o $(EXE) $(LDFLAGS)

%.o:%.cpp *.h
	$(CXX) $(CXXFLAGS) -I./cub -c $< -o $(notdir $@)

clean:
	rm -f *.o $(EXE)

realclean: clean

distclean: clean


