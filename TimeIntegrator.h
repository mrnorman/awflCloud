
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "Tendencies.h"

class TimeIntegrator {

  real4d stateTmp;
  real4d tend;
  real4d tendTmp;
  Tendencies tendencies;
  int dsSwitch;

public :


  void initialize(Domain &dom);


  void stepForward(real4d &state, Domain &dom, Exchange &exch, Parallel const &par);


  void stepForwardADER(real4d &state, Domain &dom, Exchange &exch, Parallel const &par);


  void stepForwardSSPRK3(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par);


  void applyTendencies(real4d &state2, real const c0, real4d const &state0,
                                              real const c1, real4d const &state1,
                                              real const ct, real4d const &tend, Domain const &dom);


  void appendTendencies(real4d &tend, real4d const &tendTmp, Domain const &dom);


  void applyHeatingCooling(real4d &state, Parallel const &par, Domain const &dom);

};

#endif
