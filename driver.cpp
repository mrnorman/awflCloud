
#include "stdlib.h"
#include <iostream>
#include <string>
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "State.h"
#include "Initializer.h"
#include "TimeIntegrator.h"
#include "FileIO.h"
#include "Exchange.h"
#include <Kokkos_Core.hpp>

int main(int argc, char** argv) {
  Kokkos::initialize();

  {
    // Create the model objects
    State          state;
    Domain         dom;
    Parallel       par;
    Parser         parser;
    Initializer    init;
    FileIO         io;
    Exchange       exch;
    TimeIntegrator tint;


    // Initialize MPI and read the input file
    init.initialize_mpi( &argc , &argv , par );

    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    parser.readParamsFile(inFile, dom, par, io);

    // Initialize the model
    init.initialize(state, dom, par, exch, tint);

    // Output the initial model state
    io.outputInit(state, dom, par);

    while (dom.etime < dom.simLength) {
      if (dom.etime + dom.dt > dom.simLength) { dom.dt = dom.simLength - dom.etime; }
      tint.stepForward(state, dom, exch, par);
      dom.etime += dom.dt;
      if (par.masterproc) {std::cout << dom.etime << "\n";}
      io.output(state, dom, par);
    }

    io.output(state, dom, par);
  }

  Kokkos::finalize();

}
