
#include "stdlib.h"
#include <iostream>
#include <string>
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "State.h"
#include "Initializer.h"
#include "Array.h"

int main(int argc, char** argv) {
  State state;
  Domain dom;
  Parallel par;
  Parser parser;
  Initializer init;

  // Default input file is "input.txt" unless the user passes in another file
  std::string inFile = "input.txt";
  if (argc > 1) inFile = argv[1];
  parser.readParamsFile(inFile, dom, par);

  // Initialize the model
  init.initialize(state, dom, par);
}
