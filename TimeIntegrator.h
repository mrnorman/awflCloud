
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
  int dsSwitch;

public :


  inline void initialize(Domain &dom) {
    if (timeMethod == TIME_SSPRK3) {
      stateTmp  .setup(numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
      tendArrTmp.setup(numState,dom.nz,dom.ny,dom.nx);
    }
    tendArr.setup(numState,dom.nz,dom.ny,dom.nx);
    tend.initialize(dom);
    dsSwitch = 1;
  }


  inline void stepForward(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    if (timeMethod == TIME_SSPRK3) {
      stepForwardSSPRK3(state, dom, exch, par);
    } else if (timeMethod == TIME_ADER) {
      stepForwardADER(state, dom, exch, par);
    } else {
      std::cout << "Error: Unrecognized timeMethod\n";
      exit(-1);
    }
  }


  inline void stepForwardADER(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    if (dsSwitch) {
      dsSwitch = 0;
      tend.compEulerTendADER_X(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
      if (!dom.run2d) {
        tend.compEulerTendADER_Y(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr);
        applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
      }
      tend.compEulerTendADER_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 0.5_fp , tendArr, dom);
      tend.compEulerTendADER_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 0.5_fp , tendArr, dom);
    } else {
      dsSwitch = 1;
      tend.compEulerTendADER_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 0.5_fp , tendArr, dom);
      tend.compEulerTendADER_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 0.5_fp , tendArr, dom);
      if (!dom.run2d) {
        tend.compEulerTendADER_Y(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr);
        applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
      }
      tend.compEulerTendADER_X(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
    }
  }


  inline void stepForwardSSPRK3(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    // Stage 1
    tend.compEulerTendSD_X(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr   );
    if (!dom.run2d) {
      tend.compEulerTendSD_Y(state.state, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    }
    tend.compEulerTendSD_Z(state.state, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(state.state,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( stateTmp , 1._fp , state.state , 0._fp , stateTmp , 1._fp , tendArr, dom);

    // Stage 2
    tend.compEulerTendSD_X(stateTmp, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr   );
    if (!dom.run2d) {
      tend.compEulerTendSD_Y(stateTmp, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    }
    tend.compEulerTendSD_Z(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(stateTmp,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( stateTmp , 0.75_fp , state.state , 0.25_fp , stateTmp , 0.25_fp , tendArr, dom);

    // Stage 3
    tend.compEulerTendSD_X(stateTmp, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArr   );
    if (!dom.run2d) {
      tend.compEulerTendSD_Y(stateTmp, state.hyDensCells, state.hyDensThetaCells, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    }
    tend.compEulerTendSD_Z(stateTmp, state.hyDensGLL, state.hyDensThetaGLL, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    tend.compEulerTendSD_S(stateTmp,                                        dom,            tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( state.state , 1._fp/3._fp , state.state , 2._fp/3._fp , stateTmp , 2._fp/3._fp , tendArr , dom);
  }


  inline void applyTendencies(Array<real> &state2, real const c0, Array<real> &state0,
                                                   real const c1, Array<real> &state1,
                                                   real const ct, Array<real> &tend, Domain &dom) {
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            state2(l,hs+k,hs+j,hs+i) = c0 * state0(l,hs+k,hs+j,hs+i) + c1 * state1(l,hs+k,hs+j,hs+i) + ct * dom.dt * tend(l,k,j,i);
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
