
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Array.h"
#include "Domain.h"
#include "State.h"

class TimeIntegrator {

  Array<real> stateTmp;
  Tendencies tend;

public :

  inline void initialize(Domain &dom) {
    stateTmp.setup(numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
  }

  inline void stepForward(State &state, Domain &dom, Parallel &par) {
    stepForwardSSPRK3(state, dom, par);
  }

  inline void stepForwardSSPRK3(State &state, Domain &dom, Parallel &par) {

  }

};

#endif
