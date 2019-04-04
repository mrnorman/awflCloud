
#ifndef _PARSER_H_
#define _PARSER_H_

#include <fstream>
#include <string>
#include "string.h"
#include <sstream>
#include <algorithm>
#include "Domain.h"
#include "Parallel.h"


class Parser {

public:

  void readParamsFile(std::string fNameIn, Domain &dom, Parallel &par) {

    // Initialize all read-in values to -999
    dom.nx_glob   = -999;
    dom.ny_glob   = -999;
    dom.nz_glob   = -999;
    dom.xlen      = -999;
    dom.ylen      = -999;
    dom.zlen      = -999;
    dom.cfl       = -999;
    dom.simLength = -999;
    par.px        = -999;
    par.py        = -999;

    // Read in colon-separated key: value file line by line
    std::ifstream fInStream(fNameIn);
    std::string line;
    while (std::getline(fInStream, line)) {
      // Remove spaces and tabs from the line
      line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
      line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

      if (!line.empty()) {
        uint splitloc = line.find(':',0);
        std::string key   = line.substr(0,splitloc);
        std::string value = line.substr(splitloc+1,line.length()-splitloc);

        std::stringstream ssVal(value);

        if ( !strcmp( "nx"        , key.c_str() ) ) { ssVal >> dom.nx_glob  ; }
        if ( !strcmp( "ny"        , key.c_str() ) ) { ssVal >> dom.ny_glob  ; }
        if ( !strcmp( "nz"        , key.c_str() ) ) { ssVal >> dom.nz_glob  ; }
        if ( !strcmp( "xlen"      , key.c_str() ) ) { ssVal >> dom.xlen     ; }
        if ( !strcmp( "ylen"      , key.c_str() ) ) { ssVal >> dom.ylen     ; }
        if ( !strcmp( "zlen"      , key.c_str() ) ) { ssVal >> dom.zlen     ; }
        if ( !strcmp( "cfl"       , key.c_str() ) ) { ssVal >> dom.cfl      ; }
        if ( !strcmp( "simLength" , key.c_str() ) ) { ssVal >> dom.simLength; }
        if ( !strcmp( "parNx"     , key.c_str() ) ) { ssVal >> par.px       ; }
        if ( !strcmp( "parNy"     , key.c_str() ) ) { ssVal >> par.py       ; }
      }
    }

    // Test to make sure all values were initialized
    if (dom.nx_glob   == -999) { std::cout << "Error: key " << "nx"        << " not set."; exit(-1); }
    if (dom.ny_glob   == -999) { std::cout << "Error: key " << "ny"        << " not set."; exit(-1); }
    if (dom.nz_glob   == -999) { std::cout << "Error: key " << "nz"        << " not set."; exit(-1); }
    if (dom.xlen      == -999) { std::cout << "Error: key " << "xlen"      << " not set."; exit(-1); }
    if (dom.ylen      == -999) { std::cout << "Error: key " << "ylen"      << " not set."; exit(-1); }
    if (dom.zlen      == -999) { std::cout << "Error: key " << "zlen"      << " not set."; exit(-1); }
    if (dom.cfl       == -999) { std::cout << "Error: key " << "cfl"       << " not set."; exit(-1); }
    if (dom.simLength == -999) { std::cout << "Error: key " << "simLength" << " not set."; exit(-1); }
    if (par.px        == -999) { std::cout << "Error: key " << "parNx"     << " not set."; exit(-1); }
    if (par.py        == -999) { std::cout << "Error: key " << "parNy"     << " not set."; exit(-1); }

    // Print out the values
    std::cout << "nx: "        << dom.nx_glob   << "\n";
    std::cout << "ny: "        << dom.ny_glob   << "\n";
    std::cout << "nz: "        << dom.nz_glob   << "\n";
    std::cout << "xlen: "      << dom.xlen      << "\n";
    std::cout << "ylen: "      << dom.ylen      << "\n";
    std::cout << "zlen: "      << dom.zlen      << "\n";
    std::cout << "cfl: "       << dom.cfl       << "\n";
    std::cout << "simLength: " << dom.simLength << "\n";
    std::cout << "parNx: "     << par.px        << "\n";
    std::cout << "parNy: "     << par.py        << "\n";

  }

};

#endif
