
#include "Tendencies.h"



void Tendencies::initialize(Domain const &dom) {
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


void Tendencies::compEulerTend_X(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
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


void Tendencies::compEulerTend_Y(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
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


void Tendencies::compEulerTend_Z(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
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


void Tendencies::compEulerTend_S(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
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


void Tendencies::compStrakaTend(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
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



