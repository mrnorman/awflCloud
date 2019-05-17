
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

  double dtWallTime = 0;
  double ioWallTime = 0;
  double initWallTime = 0;
  double timeTmp;

  Kokkos::Timer timer;

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

    timeTmp = timer.seconds();

    // Initialize MPI and read the input file
    init.initialize_mpi( &argc , &argv , par );

    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    parser.readParamsFile(inFile, dom, par, io);

    // Initialize the model
    init.initialize(state, dom, par, exch, tint);

    initWallTime += timer.seconds() - timeTmp;

    // Output the initial model state
    timeTmp = timer.seconds();
    io.outputInit(state, dom, par);
    ioWallTime += timer.seconds() - timeTmp;

    while (dom.etime < dom.simLength) {
      if (dom.etime + dom.dt > dom.simLength) { dom.dt = dom.simLength - dom.etime; }

      timeTmp = timer.seconds();
      tint.stepForward(state, dom, exch, par);
      dtWallTime += timer.seconds() - timeTmp;

      dom.etime += dom.dt;
      if (par.masterproc) {std::cout << dom.etime << "\n";}

      timeTmp = timer.seconds();
      io.output(state, dom, par);
      ioWallTime += timer.seconds() - timeTmp;
    }

    std::cout << "Initialization walltime: " << initWallTime << " seconds.\n";
    std::cout << "Time stepping walltime:  " << dtWallTime   << " seconds.\n";
    std::cout << "File IO walltime:        " << ioWallTime   << " seconds.\n";
  }

  Kokkos::finalize();

}
