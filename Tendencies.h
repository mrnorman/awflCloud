
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
#include "TendenciesThetaPrimADER.h"
#include "TendenciesThetaPrimSD.h"

class Tendencies {
  TendenciesThetaConsADER  tendenciesThetaConsADER;
  TendenciesThetaConsSD    tendenciesThetaConsSD  ;
  TendenciesThetaPrimADER  tendenciesThetaPrimADER;
  TendenciesThetaPrimSD    tendenciesThetaPrimSD  ;

  int const useTendThetaConsADER = 1;
  int const useTendThetaConsSD   = 2;
  int const useTendThetaPrimADER = 3;
  int const useTendThetaPrimSD   = 4;

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
      if (timeMethod == TIME_ADER) {
        tendenciesThetaPrimADER.initialize(dom);
        useTend = useTendThetaPrimADER;
      } else if (timeMethod == TIME_SSPRK3) {
        tendenciesThetaPrimSD.initialize(dom);
        useTend = useTendThetaPrimSD;
      }
    }
  }


  inline void compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_X(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_X(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimADER) {
      tendenciesThetaPrimADER.compEulerTend_X(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_X(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_Y(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_Y(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimADER) {
      tendenciesThetaPrimADER.compEulerTend_Y(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_Y(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_Z(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_Z(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimADER) {
      tendenciesThetaPrimADER.compEulerTend_Z(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_Z(state, dom, exch, par, tend);
    }
  }


  inline void compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compEulerTend_S(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compEulerTend_S(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimADER) {
      tendenciesThetaPrimADER.compEulerTend_S(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compEulerTend_S(state, dom, exch, par, tend);
    }
  }


  inline void compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend) {
    if        (useTend == useTendThetaConsADER) {
      tendenciesThetaConsADER.compStrakaTend(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaConsSD  ) {
      tendenciesThetaConsSD  .compStrakaTend(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimADER) {
      tendenciesThetaPrimADER.compStrakaTend(state, dom, exch, par, tend);
    } else if (useTend == useTendThetaPrimSD  ) {
      tendenciesThetaPrimSD  .compStrakaTend(state, dom, exch, par, tend);
    }
  }

};

#endif

