
#include "stdlib.h"
#include "iostream"
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"

int main(int argc, char** argv) {
  Domain dom;
  Parallel par;
  Parser parser;
  parser.readParamsFile("input.txt", dom, par);
}
