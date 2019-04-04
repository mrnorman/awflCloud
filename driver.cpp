
#include "stdlib.h"
#include <iostream>
#include <string>
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"

int main(int argc, char** argv) {
  Domain dom;
  Parallel par;
  Parser parser;

  std::string inFile = "input.txt";
  if (argc > 1) inFile = argv[1];
  parser.readParamsFile(inFile, dom, par);
}
