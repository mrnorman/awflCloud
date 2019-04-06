
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
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            state.state(l,k,j,hs+i) += dom.dt * tend.tend(l,k,j,i);
          }
        }
      }
    }
  }

};

#endif
