
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
    tend.compEulerTendSD_X(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par);
    applyTendencies(state.state, dom, tend.tend);

    tend.compEulerTendSD_Y(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par);
    applyTendencies(state.state, dom, tend.tend);

    tend.compEulerTendSD_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par);
    applyTendencies(state.state, dom, tend.tend);
  }

  inline void applyTendencies(Array<real> &state, Domain const &dom, Array<real> const &tend) {
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            state(l,hs+k,hs+j,hs+i) += dom.dt * tend(l,k,j,i);
          }
        }
      }
    }
  }

};

#endif
