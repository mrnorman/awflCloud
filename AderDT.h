
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"

class AderDT {


public:


  inline void diffTransformEulerX( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux, SArray<real,tord,tord> &deriv ) {
    SArray<real,tord,tord> ruu, ruv, ruw, rut, rtgamma;
    SArray<real,tord> tot_ruu, tot_ruv, tot_ruw, tot_rut, tot_rtgamma;

    // Zero out intermediate arrays
    for (int kt=1; kt<tord; kt++) {
      for (int ii=0; ii<tord; ii++) {
        ruu    (kt,ii) = 0;
        ruv    (kt,ii) = 0;
        ruw    (kt,ii) = 0;
        rut    (kt,ii) = 0;
        rtgamma(kt,ii) = 0;
      }
    }

    // Compute the zeroth-order DTs of the intermediate functions and fluxes
    for (int ii=0; ii<tord; ii++) {
      real r = state(idR ,0,ii);
      real u = state(idRU,0,ii) / r;
      real v = state(idRV,0,ii) / r;
      real w = state(idRW,0,ii) / r;
      real t = state(idTH,0,ii) / r;

      ruu    (0,ii) = r*u*u;
      ruv    (0,ii) = r*u*v;
      ruw    (0,ii) = r*u*w;
      rut    (0,ii) = r*u*t;
      rtgamma(0,ii) = mypow( r*t , GAMMA );

      flux(idR ,0,ii) = r*u;
      flux(idRU,0,ii) = ruu(0,ii) + C0*rtgamma(0,ii);
      flux(idRV,0,ii) = ruv(0,ii);
      flux(idRW,0,ii) = ruw(0,ii);
      flux(idTH,0,ii) = rut(0,ii);
    }

    // Loop over the time derivatives
    for (int kt=0; kt<tord-1; kt++) {
      // Compute the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dx = 0;
          for (int s=0; s<tord; s++) {
            d_dx += deriv(s,ii) * flux(l,kt,s);
          }
          state(l,kt+1,ii) = -d_dx/(kt+1);
        }
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_ruu(ii) = 0;
        tot_ruv(ii) = 0;
        tot_ruw(ii) = 0;
        tot_rut(ii) = 0;
      }
      for (int rt=0; rt<kt+1; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_ruu(ii) += state(idRU,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * ruu(kt+1-rt,ii);
          tot_ruv(ii) += state(idRU,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * ruv(kt+1-rt,ii);
          tot_ruw(ii) += state(idRU,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * ruw(kt+1-rt,ii);
          tot_rut(ii) += state(idRU,rt,ii) * state(idTH,kt+1-rt,ii) - state(idR,rt,ii) * rut(kt+1-rt,ii);
        }
      }
      for (int ii=0; ii<tord; ii++) {
        ruu(kt+1,ii) = tot_ruu(ii) / state(idR,0,ii);
        ruv(kt+1,ii) = tot_ruv(ii) / state(idR,0,ii);
        ruw(kt+1,ii) = tot_ruw(ii) / state(idR,0,ii);
        rut(kt+1,ii) = tot_rut(ii) / state(idR,0,ii);
      }

      // Compute rtgamma at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_rtgamma(ii) = 0;
      }
      for (int rt=0; rt<kt; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_rtgamma(ii) += (kt+1-rt) * ( GAMMA*rtgamma(rt,ii)*state(idTH,kt+1-rt,ii) - state(idTH,rt,ii)*rtgamma(kt+1-rt,ii) );
        }
      }
      for (int ii=0; ii<tord; ii++) {
        rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idTH,kt+1,ii) + tot_rtgamma(ii) / (kt+1) ) / state(idTH,0,ii);
      }

      // Compute the fluxes at the next time level
      for (int ii=0; ii<tord; ii++) {
        flux(idR ,kt+1,ii) = state(idRU,kt+1,ii);
        flux(idRU,kt+1,ii) = ruu(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
        flux(idRV,kt+1,ii) = ruv(kt+1,ii);
        flux(idRW,kt+1,ii) = ruw(kt+1,ii);
        flux(idTH,kt+1,ii) = rut(kt+1,ii);
      }
    }
  }


  inline void diffTransformEulerY( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux, SArray<real,tord,tord> &deriv ) {
    SArray<real,tord,tord> rvu, rvv, rvw, rvt, rtgamma;
    SArray<real,tord> tot_rvu, tot_rvv, tot_rvw, tot_rvt, tot_rtgamma;

    // Zero out intermediate arrays
    for (int kt=1; kt<tord; kt++) {
      for (int ii=0; ii<tord; ii++) {
        rvu    (kt,ii) = 0;
        rvv    (kt,ii) = 0;
        rvw    (kt,ii) = 0;
        rvt    (kt,ii) = 0;
        rtgamma(kt,ii) = 0;
      }
    }

    // Compute the zeroth-order DTs of the intermediate functions and fluxes
    for (int ii=0; ii<tord; ii++) {
      real r = state(idR ,0,ii);
      real u = state(idRU,0,ii) / r;
      real v = state(idRV,0,ii) / r;
      real w = state(idRW,0,ii) / r;
      real t = state(idTH,0,ii) / r;

      rvu    (0,ii) = r*v*u;
      rvv    (0,ii) = r*v*v;
      rvw    (0,ii) = r*v*w;
      rvt    (0,ii) = r*v*t;
      rtgamma(0,ii) = mypow( r*t , GAMMA );

      flux(idR ,0,ii) = r*v;
      flux(idRU,0,ii) = rvu(0,ii);
      flux(idRV,0,ii) = rvv(0,ii) + C0*rtgamma(0,ii);
      flux(idRW,0,ii) = rvw(0,ii);
      flux(idTH,0,ii) = rvt(0,ii);
    }

    // Loop over the time derivatives
    for (int kt=0; kt<tord-1; kt++) {
      // Compute the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dx = 0;
          for (int s=0; s<tord; s++) {
            d_dx += deriv(s,ii) * flux(l,kt,s);
          }
          state(l,kt+1,ii) = -d_dx/(kt+1);
        }
      }

      // Compute rv* at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_rvu(ii) = 0;
        tot_rvv(ii) = 0;
        tot_rvw(ii) = 0;
        tot_rvt(ii) = 0;
      }
      for (int rt=0; rt<kt+1; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_rvu(ii) += state(idRV,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * rvu(kt+1-rt,ii);
          tot_rvv(ii) += state(idRV,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * rvv(kt+1-rt,ii);
          tot_rvw(ii) += state(idRV,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * rvw(kt+1-rt,ii);
          tot_rvt(ii) += state(idRV,rt,ii) * state(idTH,kt+1-rt,ii) - state(idR,rt,ii) * rvt(kt+1-rt,ii);
        }
      }
      for (int ii=0; ii<tord; ii++) {
        rvu(kt+1,ii) = tot_rvu(ii) / state(idR,0,ii);
        rvv(kt+1,ii) = tot_rvv(ii) / state(idR,0,ii);
        rvw(kt+1,ii) = tot_rvw(ii) / state(idR,0,ii);
        rvt(kt+1,ii) = tot_rvt(ii) / state(idR,0,ii);
      }

      // Compute rtgamma at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_rtgamma(ii) = 0;
      }
      for (int rt=0; rt<kt; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_rtgamma(ii) += (kt+1-rt) * ( GAMMA*rtgamma(rt,ii)*state(idTH,kt+1-rt,ii) - state(idTH,rt,ii)*rtgamma(kt+1-rt,ii) );
        }
      }
      for (int ii=0; ii<tord; ii++) {
        rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idTH,kt+1,ii) + tot_rtgamma(ii) / (kt+1) ) / state(idTH,0,ii);
      }

      // Compute the fluxes at the next time level
      for (int ii=0; ii<tord; ii++) {
        flux(idR ,kt+1,ii) = state(idRV,kt+1,ii);
        flux(idRU,kt+1,ii) = rvu(kt+1,ii);
        flux(idRV,kt+1,ii) = rvv(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
        flux(idRW,kt+1,ii) = rvw(kt+1,ii);
        flux(idTH,kt+1,ii) = rvt(kt+1,ii);
      }
    }
  }


  inline void diffTransformEulerZ( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux,
                                   SArray<real,numState,tord,tord> &source, SArray<real,tord,tord> &deriv, SArray<real,tord> &hyRHOT, SArray<real,tord> &hyRHO ) {
    SArray<real,tord,tord> rwu, rwv, rww, rwt, rtgamma;
    SArray<real,tord> tot_rwu, tot_rwv, tot_rww, tot_rwt, tot_rtgamma;

    // Zero out intermediate arrays
    for (int kt=1; kt<tord; kt++) {
      for (int ii=0; ii<tord; ii++) {
        rwu    (kt,ii) = 0;
        rwv    (kt,ii) = 0;
        rww    (kt,ii) = 0;
        rwt    (kt,ii) = 0;
        rtgamma(kt,ii) = 0;
      }
    }
    source = 0;

    // Compute the zeroth-order DTs of the intermediate functions and fluxes
    for (int ii=0; ii<tord; ii++) {
      real r = state(idR ,0,ii);
      real u = state(idRU,0,ii) / r;
      real v = state(idRV,0,ii) / r;
      real w = state(idRW,0,ii) / r;
      real t = state(idTH,0,ii) / r;

      rwu    (0,ii) = r*w*u;
      rwv    (0,ii) = r*w*v;
      rww    (0,ii) = r*w*w;
      rwt    (0,ii) = r*w*t;
      rtgamma(0,ii) = mypow( r*t , GAMMA );

      flux(idR ,0,ii) = r*w;
      flux(idRU,0,ii) = rwu(0,ii);
      flux(idRV,0,ii) = rwv(0,ii);
      flux(idRW,0,ii) = rww(0,ii) + C0*rtgamma(0,ii) - C0*mypow(hyRHOT(ii),GAMMA);
      flux(idTH,0,ii) = rwt(0,ii);

      source(idRW,0,ii) = -( state(idR,0,ii) - hyRHO(ii) )*GRAV;
    }

    // Loop over the time derivatives
    for (int kt=0; kt<tord-1; kt++) {
      // Compute the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dx = 0;
          for (int s=0; s<tord; s++) {
            d_dx += deriv(s,ii) * flux(l,kt,s);
          }
          state(l,kt+1,ii) = ( -d_dx + source(l,kt,ii) ) / (kt+1);
        }
      }

      // Compute rw* at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_rwu(ii) = 0;
        tot_rwv(ii) = 0;
        tot_rww(ii) = 0;
        tot_rwt(ii) = 0;
      }
      for (int rt=0; rt<kt+1; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_rwu(ii) += state(idRW,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * rwu(kt+1-rt,ii);
          tot_rwv(ii) += state(idRW,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * rwv(kt+1-rt,ii);
          tot_rww(ii) += state(idRW,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * rww(kt+1-rt,ii);
          tot_rwt(ii) += state(idRW,rt,ii) * state(idTH,kt+1-rt,ii) - state(idR,rt,ii) * rwt(kt+1-rt,ii);
        }
      }
      for (int ii=0; ii<tord; ii++) {
        rwu(kt+1,ii) = tot_rwu(ii) / state(idR,0,ii);
        rwv(kt+1,ii) = tot_rwv(ii) / state(idR,0,ii);
        rww(kt+1,ii) = tot_rww(ii) / state(idR,0,ii);
        rwt(kt+1,ii) = tot_rwt(ii) / state(idR,0,ii);
      }

      // Compute rtgamma at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_rtgamma(ii) = 0;
      }
      for (int rt=0; rt<kt; rt++) {
        for (int ii=0; ii<tord; ii++) {
          tot_rtgamma(ii) += (kt+1-rt) * ( GAMMA*rtgamma(rt,ii)*state(idTH,kt+1-rt,ii) - state(idTH,rt,ii)*rtgamma(kt+1-rt,ii) );
        }
      }
      for (int ii=0; ii<tord; ii++) {
        rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idTH,kt+1,ii) + tot_rtgamma(ii) / (kt+1) ) / state(idTH,0,ii);
      }

      // Compute the fluxes at the next time level
      for (int ii=0; ii<tord; ii++) {
        flux(idR ,kt+1,ii) = state(idRW,kt+1,ii);
        flux(idRU,kt+1,ii) = rwu(kt+1,ii);
        flux(idRV,kt+1,ii) = rwv(kt+1,ii);
        flux(idRW,kt+1,ii) = rww(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
        flux(idTH,kt+1,ii) = rwt(kt+1,ii);
        source(idRW,kt+1,ii) = -state(idR,kt+1,ii)*GRAV;
      }
    }
  }


  inline void timeAvg( SArray<real,numState,tord,tord> &dts , Domain &dom ) {
    real dtmult = dom.dt;
    for (int kt=1; kt<tord; kt++) {
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1);
        }
      }
      dtmult *= dom.dt;
    }
  }

};

#endif
