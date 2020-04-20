
#ifndef _TENDENCIES_THETA_CONS_SD_H_
#define _TENDENCIES_THETA_CONS_SD_H_

#include "const.h"
#include "Parallel.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"
#include "WenoLimiter.h"
#include "AderDT.h"
#include "TransformMatrices.h"

class TendenciesThetaConsSD {

  real5d stateLimits;
  real5d fluxLimits;
  real4d flux;
  real3d src;
  real5d stateGLL;
  SArray<real,1,tord> gllWts;
  SArray<real,2,ord,tord> to_gll;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> wenoIdl;
  real wenoSigma;

public :


  void initialize(Domain const &dom);


  // Transform ord stencil cell averages into tord GLL point values
  YAKL_INLINE void reconStencil(SArray<real,1,ord> const &stencil, SArray<real,1,tord> &gll, int const doWeno,
                                SArray<real,3,ord,ord,ord> const &wenoRecon, SArray<real,2,ord,tord> const &to_gll,
                                SArray<real,1,hs+2> const &wenoIdl, real wenoSigma) {
    SArray<real,1,ord> coefs;
    if (doWeno) {
      compute_weno_coefs(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);
    } else {
      for (int ii=0; ii<ord; ii++) {
        coefs(ii) = stencil(ii);
      }
    }

    for (int ii=0; ii<tord; ii++) {
      gll(ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(ii) += to_gll(s,ii) * coefs(s);
      }
    }
  }


  void compEulerTend_X(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend);

  void compEulerTend_Y(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend);

  void compEulerTend_Z(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend);

  void compEulerTend_S(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend);

  void stateBoundariesZ(real4d &state, Domain const &dom);

  void edgeBoundariesZ(real5d &stateLimits, real5d &fluxLimits, Domain const &dom);

  void compStrakaTend(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend);

};

#endif

