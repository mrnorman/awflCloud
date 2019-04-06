
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Array.h"
#include "Domain.h"
#include "State.h"
#include "Tendencies.h"

class TimeIntegrator {

  Array<real> stateTmp;
  Tendencies tend;

public :

  inline void initialize(Domain &dom) {
    stateTmp.setup(numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
    tend.initialize(dom);
  }

  inline void stepForward(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    stepForwardSSPRK3(state, dom, exch, par);
  }

  inline void stepForwardSSPRK3(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    tend.compEulerTendSD_X(state.state, dom, exch, par);
  }

};

#endif
