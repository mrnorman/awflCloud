
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
  Array<real> tendArr;
  Array<real> tendArrTmp;
  Tendencies tend;

public :

  inline void initialize(Domain &dom) {
    stateTmp  .setup(numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
    tendArr   .setup(numState,dom.nz,dom.ny,dom.nx);
    tendArrTmp.setup(numState,dom.nz,dom.ny,dom.nx);
    tend.initialize(dom);
  }

  inline void stepForward(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    stepForwardSSPRK3(state, dom, exch, par);
  }

  inline void stepForwardSSPRK3(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    // Stage 1
    tend.compEulerTendSD_X(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr   );
    tend.compEulerTendSD_Y(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(state.state,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies(stateTmp, state.state, dom, tendArr, dom.dt/3);

    // Stage 2
    tend.compEulerTendSD_X(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr   );
    tend.compEulerTendSD_Y(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_Z(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(stateTmp,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies(stateTmp, state.state, dom, tendArr, dom.dt/2);

    // Stage 3
    tend.compEulerTendSD_X(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr   );
    tend.compEulerTendSD_Y(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_Z(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(stateTmp,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies(state.state, state.state, dom, tendArr, dom.dt/1);
  }

  inline void applyTendencies(Array<real> &state1, Array<real> &state0, Domain const &dom, Array<real> const &tend, real const dt) {
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            state1(l,hs+k,hs+j,hs+i) = state0(l,hs+k,hs+j,hs+i) + dt * tend(l,k,j,i);
          }
        }
      }
    }
  }

  inline void appendTendencies(Array<real> &tend, Array<real> &tendTmp, Domain &dom) {
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            tend(l,k,j,i) += tendTmp(l,k,j,i);
          }
        }
      }
    }
  }

};

#endif
