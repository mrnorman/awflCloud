
#ifndef _STATE_H_
#define _STATE_H_

#include "const.h"
#include "SArray.h"

class State {

public:

  // Fluid state
  real4d state;

  // Hydrostatic background state cell averages (no halos)
  real1d hyDensCells;
  real1d hyDensThetaCells;
  real1d hyPressureCells;

  // Hydrostatic background state at tord GLL points within cells (no halos)
  real2d hyDensGLL;
  real2d hyDensThetaGLL;
  real2d hyPressureGLL;

};

#endif
