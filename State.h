
#ifndef _STATE_H_
#define _STATE_H_

#include "const.h"
#include "SArray.h"
#include "Array.h"

class State {

public:

  // Fluid state
  Array<real> state;

  // Hydrostatic background state cell averages with halos
  Array<real> hyDensCells;
  Array<real> hyDensThetaCells;
  Array<real> hyPressureCells;

  // Hydrostatic background state at tord GLL points within cells (no halos)
  Array<real> hyDensGLL;
  Array<real> hyDensThetaGLL;
  Array<real> hyPressureGLL;

};

#endif
