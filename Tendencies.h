
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"
#include "TendenciesThetaConsADER.h"
#include "TendenciesThetaConsSD.h"
#include "TendenciesThetaPrimSD.h"

class Tendencies {
  TendenciesThetaConsADER  tendenciesThetaConsADER;
  TendenciesThetaConsSD    tendenciesThetaConsSD  ;
  TendenciesThetaPrimSD    tendenciesThetaPrimSD  ;

  int const useTendThetaConsADER = 1;
  int const useTendThetaConsSD   = 2;
  int const useTendThetaPrimSD   = 3;

  int useTend;

public :


  inline void initialize(Domain const &dom) {
    if (dom.eqnSet == EQN_THETA_CONS) {
      if (timeMethod == TIME_ADER) {
        tendenciesThetaConsADER.initialize(dom);
        useTend = useTendThetaConsADER;
      } else if (timeMethod == TIME_SSPRK3) {
        tendenciesThetaConsSD.initialize(dom);
        useTend = useTendThetaConsSD;
      }
    } else if (dom.eqnSet == EQN_THETA_PRIM) {
      if (timeMethod == TIME_SSPRK3) {
        tendenciesThetaPrimSD.initialize(dom);
        useTend = useTendThetaPrimSD;
      }
    }
  }


  inline void compEulerTend_X(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_X(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_X(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_X(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_Y(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_Y(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_Y(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_Y(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_Z(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_Z(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_Z(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_Z(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_S(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_S(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_S(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_S(state, dom, exch, par, tend);
    }
  }


  inline void compStrakaTend(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compStrakaTend(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compStrakaTend(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compStrakaTend(state, dom, exch, par, tend);
    }
  }

};

#endif

