
#pragma once

#include "const.h"

class Domain {

public:

  ulong nx_glob;
  ulong ny_glob;
  ulong nz_glob;

  int nx;
  int ny;
  int nz;

  int run2d;

  int doWeno;

  real xlen;
  real ylen;
  real zlen;

  real dx;
  real dy;
  real dz;

  real cfl;
  real simLength;

  real etime;

  int eqnSet;
  int dataInit;

  real dt;

  // Hydrostatic background state cell averages (no halos)
  real1d hyDensCells;
  real1d hyDensThetaCells;
  real1d hyThetaCells;
  real1d hyPressureCells;

  // Hydrostatic background state at tord GLL points within cells (no halos)
  real2d hyDensGLL;
  real2d hyDensThetaGLL;
  real2d hyThetaGLL;
  real2d hyPressureGLL;
};

