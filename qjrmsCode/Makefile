include mach.inc

EXE = cloudFV
SOURCES = driver.cpp FileIO.cpp params.cpp Exchange.cpp Initializer.cpp Parser.cpp Tendencies.cpp TimeIntegrator.cpp \
          TendenciesThetaConsADER.cpp TendenciesThetaConsSD.cpp TendenciesThetaPrimADER.cpp TendenciesThetaPrimSD.cpp \
          YAKL/YAKL.cpp YAKL/BuddyAllocator.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

default: main

main: $(EXE)

$(EXE): $(OBJECTS)
	$(LINK) *.o $(LIB) -o $(EXE) $(LDFLAGS)

%.o:%.cpp *.h
	$(CXX) $(CXXFLAGS) -I./cub -I./YAKL -c $< -o $(notdir $@)

clean:
	rm -f *.o $(EXE)

realclean: clean

distclean: clean


