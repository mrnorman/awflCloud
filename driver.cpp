
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
#include "Array.h"
#include "Exchange.h"

int main(int argc, char** argv) {
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
  parser.readParamsFile(inFile, dom, par);

  // Initialize the model
  init.initialize(state, dom, par, exch, tint);

  // Output the initial model state
  io.outputInit(state, dom, par);

  tint.stepForward(state, dom, exch, par);

  io.output(state, dom, par);

}
