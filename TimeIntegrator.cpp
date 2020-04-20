
#include "TimeIntegrator.h"


void TimeIntegrator::initialize(Domain &dom) {
  if (timeMethod == TIME_SSPRK3) {
    stateTmp = real4d("stateTmp",numState,dom.nz+2*hs,dom.ny+2*hs,dom.nx+2*hs);
    tendTmp  = real4d("tendTmp" ,numState,dom.nz,dom.ny,dom.nx);
  }
  tend = real4d("tend",numState,dom.nz,dom.ny,dom.nx);
  tendencies.initialize(dom);
  dsSwitch = 1;
}


void TimeIntegrator::stepForward(real4d &state, Domain &dom, Exchange &exch, Parallel const &par) {
  if (timeMethod == TIME_SSPRK3) {
    stepForwardSSPRK3(state, dom, exch, par);
  } else if (timeMethod == TIME_ADER) {
    stepForwardADER  (state, dom, exch, par);
  } else {
    std::cout << "Error: Unrecognized timeMethod\n";
    exit(-1);
  }
  if (strakaVisc) {
    tendencies.compStrakaTend(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
  }
}


void TimeIntegrator::stepForwardADER(real4d &state, Domain &dom, Exchange &exch, Parallel const &par) {
  if (dsSwitch) {
    dsSwitch = 0;
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    }
    dom.dt /= 2;
    tendencies.compEulerTend_Z(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    tendencies.compEulerTend_Z(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    dom.dt *= 2;
  } else {
    dsSwitch = 1;
    dom.dt /= 2;
    tendencies.compEulerTend_Z(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    tendencies.compEulerTend_Z(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    dom.dt *= 2;
    if (!dom.run2d) {
      tendencies.compEulerTend_Y(state, dom, exch, par, tend);
      applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
    }
    tendencies.compEulerTend_X(state, dom, exch, par, tend);
    applyTendencies( state , 1._fp , state , 0._fp , state , 1._fp , tend, dom);
  }
  // applyHeatingCooling(state,par,dom);
}


void TimeIntegrator::stepForwardSSPRK3(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par) {
  // Stage 1
  tendencies.compEulerTend_X(state, dom, exch, par, tend   );
  if (!dom.run2d) {
    tendencies.compEulerTend_Y(state, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  }
  tendencies.compEulerTend_Z(state, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  tendencies.compEulerTend_S(state, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  applyTendencies( stateTmp , 1._fp , state , 0._fp , stateTmp , 1._fp , tend, dom);

  // Stage 2
  tendencies.compEulerTend_X(stateTmp, dom, exch, par, tend   );
  if (!dom.run2d) {
    tendencies.compEulerTend_Y(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  }
  tendencies.compEulerTend_Z(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  tendencies.compEulerTend_S(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  applyTendencies( stateTmp , 0.75_fp , state , 0.25_fp , stateTmp , 0.25_fp , tend, dom);

  // Stage 3
  tendencies.compEulerTend_X(stateTmp, dom, exch, par, tend   );
  if (!dom.run2d) {
    tendencies.compEulerTend_Y(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  }
  tendencies.compEulerTend_Z(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  tendencies.compEulerTend_S(stateTmp, dom, exch, par, tendTmp); appendTendencies(tend, tendTmp, dom);
  applyTendencies( state , 1._fp/3._fp , state , 2._fp/3._fp , stateTmp , 2._fp/3._fp , tend , dom);
}


void TimeIntegrator::applyTendencies(real4d &state2, real const c0, real4d const &state0,
                                            real const c1, real4d const &state1,
                                            real const ct, real4d const &tend, Domain const &dom) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  parallel_for( Bounds<4>(numState,dom.nz,dom.ny,dom.nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
    state2(l,hs+k,hs+j,hs+i) = c0 * state0(l,hs+k,hs+j,hs+i) + c1 * state1(l,hs+k,hs+j,hs+i) + ct * dom.dt * tend(l,k,j,i);
  });
}


void TimeIntegrator::appendTendencies(real4d &tend, real4d const &tendTmp, Domain const &dom) {
  // for (int l=0; l<numState; l++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int i=0; i<dom.nx; i++) {
  parallel_for( Bounds<4>(numState,dom.nz,dom.ny,dom.nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
    tend(l,k,j,i) += tendTmp(l,k,j,i);
  });
}


void TimeIntegrator::applyHeatingCooling(real4d &state, Parallel const &par, Domain const &dom) {
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  parallel_for( Bounds<2>(dom.ny,dom.nx) , YAKL_LAMBDA (int j, int i) {
    real xloc = (par.i_beg + i + 0.5_fp) * dom.dx;
    real yloc = (par.j_beg + j + 0.5_fp) * dom.dy;
    if (dom.run2d) {yloc = dom.ylen/2;}
    for (int k=dom.nz; k < hs+dom.nz; k++) {
      state(idRT,k,hs+j,hs+i) -= dom.dt * 0.01;
    }   
    for (int k=hs; k < 2*hs; k++) {
      state(idRT,k,hs+j,hs+i) += dom.dt * 0.01;
    }   
  });
}


