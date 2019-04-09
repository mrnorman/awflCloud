
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Array.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"
#include "TransformMatrices.h"

class Tendencies {

  Array<real> stateLimits;
  Array<real> fluxLimits;
  Array<real> flux;
  TransformMatrices<real> trans;
  Riemann riem;

public :


  inline void initialize(Domain &dom) {
    fluxLimits .setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    stateLimits.setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    flux       .setup(numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
  }


  inline void compEulerTendSD_X(Array<real> &state, Array<real> &hyDensCells, Array<real> &hyDensThetaCells,
                                Domain &dom, Exchange &exch, Parallel &par, Array<real> &tend) {
    SArray<real,ord,ord,ord> s2g_lower_tmp;
    SArray<real,ord,tord> s2g_lower;

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    trans.sten_to_gll_lower( s2g_lower_tmp );
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        s2g_lower(j,i) = s2g_lower_tmp(tord-1,j,i);
      }
    }

    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState,tord> gllState;
          SArray<real,numState,tord> gllFlux;

          // Compute GLL points from cell averages
          for (int l=0; l<numState; l++) {
            for (int ii=0; ii<tord; ii++) {
              gllState(l,ii) = 0.;
              for (int s=0; s<ord; s++) {
                gllState(l,ii) += s2g_lower(s,ii) * state(l,hs+k,hs+j,i+s);
              }
            }
          }
          for (int ii=0; ii<tord; ii++) {
            gllState(idR ,ii) += hyDensCells     (hs+k);
            gllState(idTH,ii) += hyDensThetaCells(hs+k);
          }

          // Compute fluxes and at the GLL points
          for (int ii=0; ii<tord; ii++) {
            real r = gllState(idR ,ii);
            real u = gllState(idRU,ii) / r;
            real v = gllState(idRV,ii) / r;
            real w = gllState(idRW,ii) / r;
            real t = gllState(idTH,ii) / r;
            real p = C0 * mypow( r*t , GAMMA );

            gllFlux(idR ,ii) = r*u;
            gllFlux(idRU,ii) = r*u*u + p;
            gllFlux(idRV,ii) = r*u*v;
            gllFlux(idRW,ii) = r*u*w;
            gllFlux(idTH,ii) = r*u*t;
          }

          // Store state and flux limits into a globally indexed array
          for (int l=0; l<numState; l++) {
            // Store the left cell edge state and flux estimates
            stateLimits(l,1,k,j,i  ) = gllState(l,0);
            fluxLimits (l,1,k,j,i  ) = gllFlux (l,0);

            // Store the Right cell edge state and flux estimates
            stateLimits(l,0,k,j,i+1) = gllState(l,tord-1);
            fluxLimits (l,0,k,j,i+1) = gllFlux (l,tord-1);
          }

        }
      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_x   (dom, stateLimits, numState);
    exch.edgePackN_x   (dom, fluxLimits , numState);
    exch.edgeExchange_x(dom, par);
    exch.edgeUnpackN_x (dom, stateLimits, numState);
    exch.edgeUnpackN_x (dom, fluxLimits , numState);

    // Local lax-friedrichs fluxes
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx+1; i++) {
          SArray<real,numState> s1, s2, f1, f2, upw;
          for (int l=0; l<numState; l++) {
            s1(l) = stateLimits(l,0,k,j,i);
            s2(l) = stateLimits(l,1,k,j,i);
            f1(l) = fluxLimits (l,0,k,j,i);
            f2(l) = fluxLimits (l,1,k,j,i);
          }
          riem.riemannX(s1, s2, f1, f2, upw);
          for (int l=0; l<numState; l++) {
            flux(l,k,j,i) = upw(l);
          }
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            tend(l,k,j,i) = - ( flux(l,k,j,i+1) - flux(l,k,j,i) ) / dom.dx;
          }
        }
      }
    }
  }


  inline void compEulerTendSD_Y(Array<real> &state, Array<real> &hyDensCells, Array<real> &hyDensThetaCells,
                                Domain &dom, Exchange &exch, Parallel &par, Array<real> &tend) {
    SArray<real,ord,ord,ord> s2g_lower_tmp;
    SArray<real,ord,tord> s2g_lower;

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    trans.sten_to_gll_lower( s2g_lower_tmp );
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        s2g_lower(j,i) = s2g_lower_tmp(tord-1,j,i);
      }
    }

    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_y   (dom, state, numState);
    exch.haloExchange_y(dom, par);
    exch.haloUnpackN_y (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState,tord> gllState;
          SArray<real,numState,tord> gllFlux;

          // Compute GLL points from cell averages
          for (int l=0; l<numState; l++) {
            for (int ii=0; ii<tord; ii++) {
              gllState(l,ii) = 0.;
              for (int s=0; s<ord; s++) {
                gllState(l,ii) += s2g_lower(s,ii) * state(l,hs+k,j+s,hs+i);
              }
            }
          }
          for (int ii=0; ii<tord; ii++) {
            gllState(idR ,ii) += hyDensCells     (hs+k);
            gllState(idTH,ii) += hyDensThetaCells(hs+k);
          }

          // Compute fluxes and at the GLL points
          for (int ii=0; ii<tord; ii++) {
            real r = gllState(idR ,ii);
            real u = gllState(idRU,ii) / r;
            real v = gllState(idRV,ii) / r;
            real w = gllState(idRW,ii) / r;
            real t = gllState(idTH,ii) / r;
            real p = C0 * mypow( r*t , GAMMA );

            gllFlux(idR ,ii) = r*v;
            gllFlux(idRU,ii) = r*v*u;
            gllFlux(idRV,ii) = r*v*v + p;
            gllFlux(idRW,ii) = r*v*w;
            gllFlux(idTH,ii) = r*v*t;
          }

          // Store state and flux limits into a globally indexed array
          for (int l=0; l<numState; l++) {
            // Store the left cell edge state and flux estimates
            stateLimits(l,1,k,j  ,i) = gllState(l,0);
            fluxLimits (l,1,k,j  ,i) = gllFlux (l,0);

            // Store the Right cell edge state and flux estimates
            stateLimits(l,0,k,j+1,i) = gllState(l,tord-1);
            fluxLimits (l,0,k,j+1,i) = gllFlux (l,tord-1);
          }

        }
      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_y   (dom, stateLimits, numState);
    exch.edgePackN_y   (dom, fluxLimits , numState);
    exch.edgeExchange_y(dom, par);
    exch.edgeUnpackN_y (dom, stateLimits, numState);
    exch.edgeUnpackN_y (dom, fluxLimits , numState);

    // Local lax-friedrichs fluxes
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny+1; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState> s1, s2, f1, f2, upw;
          for (int l=0; l<numState; l++) {
            s1(l) = stateLimits(l,0,k,j,i);
            s2(l) = stateLimits(l,1,k,j,i);
            f1(l) = fluxLimits (l,0,k,j,i);
            f2(l) = fluxLimits (l,1,k,j,i);
          }
          riem.riemannY(s1, s2, f1, f2, upw);
          for (int l=0; l<numState; l++) {
            flux(l,k,j,i) = upw(l);
          }
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            tend(l,k,j,i) = - ( flux(l,k,j+1,i) - flux(l,k,j,i) ) / dom.dy;
          }
        }
      }
    }
  }


  inline void compEulerTendSD_Z(Array<real> &state, Array<real> &hyDensGLL, Array<real> &hyDensThetaGLL,
                                Domain &dom, Exchange &exch, Parallel &par, Array<real> &tend) {
    SArray<real,ord,ord,ord> s2g_lower_tmp;
    SArray<real,ord,tord> s2g_lower;

    // Setup the matrix to transform a stencil of ord cell averages into 2 GLL points
    trans.sten_to_gll_lower( s2g_lower_tmp );
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        s2g_lower(j,i) = s2g_lower_tmp(tord-1,j,i);
      }
    }

    // Boundaries for the fluid state in the z-direction
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        for (int ii=0; ii<hs; ii++) {
          state(idR ,ii,hs+j,hs+i) = state(idR ,hs,hs+j,hs+i);
          state(idRU,ii,hs+j,hs+i) = state(idRU,hs,hs+j,hs+i);
          state(idRV,ii,hs+j,hs+i) = state(idRV,hs,hs+j,hs+i);
          state(idRW,ii,hs+j,hs+i) = 0;
          state(idTH,ii,hs+j,hs+i) = state(idTH,hs,hs+j,hs+i);

          state(idR ,dom.nz+hs+ii,hs+j,hs+i) = state(idR ,dom.nz+hs-1,hs+j,hs+i);
          state(idRU,dom.nz+hs+ii,hs+j,hs+i) = state(idRU,dom.nz+hs-1,hs+j,hs+i);
          state(idRV,dom.nz+hs+ii,hs+j,hs+i) = state(idRV,dom.nz+hs-1,hs+j,hs+i);
          state(idRW,dom.nz+hs+ii,hs+j,hs+i) = 0;
          state(idTH,dom.nz+hs+ii,hs+j,hs+i) = state(idTH,dom.nz+hs-1,hs+j,hs+i);
        }
      }
    }

    // Reconstruct to 2 GLL points in the x-direction
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState,tord> gllState;
          SArray<real,numState,tord> gllFlux;

          // Compute GLL points from cell averages
          for (int l=0; l<numState; l++) {
            for (int ii=0; ii<tord; ii++) {
              gllState(l,ii) = 0.;
              for (int s=0; s<ord; s++) {
                gllState(l,ii) += s2g_lower(s,ii) * state(l,k+s,hs+j,hs+i);
              }
            }
          }
          for (int ii=0; ii<tord; ii++) {
            gllState(idR ,ii) += hyDensGLL     (k,ii);
            gllState(idTH,ii) += hyDensThetaGLL(k,ii);
          }

          // Boundary conditions
          if (k == 0       ) { gllState(idRW,0     ) = 0; }
          if (k == dom.nz-1) { gllState(idRW,tord-1) = 0; }

          // Compute fluxes and at the GLL points
          for (int ii=0; ii<tord; ii++) {
            real r = gllState(idR ,ii);
            real u = gllState(idRU,ii) / r;
            real v = gllState(idRV,ii) / r;
            real w = gllState(idRW,ii) / r;
            real t = gllState(idTH,ii) / r;
            real p = C0 * mypow( r*t , GAMMA );

            gllFlux(idR ,ii) = r*w;
            gllFlux(idRU,ii) = r*w*u;
            gllFlux(idRV,ii) = r*w*v;
            gllFlux(idRW,ii) = r*w*w + p - C0*mypow(hyDensThetaGLL(k,ii),GAMMA);
            gllFlux(idTH,ii) = r*w*t;
          }

          // Store state and flux limits into a globally indexed array
          for (int l=0; l<numState; l++) {
            // Store the left cell edge state and flux estimates
            stateLimits(l,1,k  ,j,i) = gllState(l,0);
            fluxLimits (l,1,k  ,j,i) = gllFlux (l,0);

            // Store the Right cell edge state and flux estimates
            stateLimits(l,0,k+1,j,i) = gllState(l,tord-1);
            fluxLimits (l,0,k+1,j,i) = gllFlux (l,tord-1);
          }

        }
      }
    }

    // Apply boundary conditions to fluxes and state values
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        stateLimits(idR ,0,0     ,j,i) = stateLimits(idR ,1,0     ,j,i);
        stateLimits(idRU,0,0     ,j,i) = stateLimits(idRU,1,0     ,j,i);
        stateLimits(idRV,0,0     ,j,i) = stateLimits(idRV,1,0     ,j,i);
        stateLimits(idRW,0,0     ,j,i) = 0;
        stateLimits(idRW,1,0     ,j,i) = 0;
        stateLimits(idTH,0,0     ,j,i) = stateLimits(idTH,1,0     ,j,i);

        stateLimits(idR ,1,dom.nz,j,i) = stateLimits(idR ,0,dom.nz,j,i);
        stateLimits(idRU,1,dom.nz,j,i) = stateLimits(idRU,0,dom.nz,j,i);
        stateLimits(idRV,1,dom.nz,j,i) = stateLimits(idRV,0,dom.nz,j,i);
        stateLimits(idRW,0,dom.nz,j,i) = 0;
        stateLimits(idRW,1,dom.nz,j,i) = 0;
        stateLimits(idTH,1,dom.nz,j,i) = stateLimits(idTH,0,dom.nz,j,i);

        fluxLimits(idR ,0,0     ,j,i) = fluxLimits(idR ,1,0     ,j,i);
        fluxLimits(idRU,0,0     ,j,i) = fluxLimits(idRU,1,0     ,j,i);
        fluxLimits(idRV,0,0     ,j,i) = fluxLimits(idRV,1,0     ,j,i);
        fluxLimits(idRW,0,0     ,j,i) = fluxLimits(idRW,1,0     ,j,i);
        fluxLimits(idTH,0,0     ,j,i) = fluxLimits(idTH,1,0     ,j,i);

        fluxLimits(idR ,1,dom.nz,j,i) = fluxLimits(idR ,0,dom.nz,j,i);
        fluxLimits(idRU,1,dom.nz,j,i) = fluxLimits(idRU,0,dom.nz,j,i);
        fluxLimits(idRV,1,dom.nz,j,i) = fluxLimits(idRV,0,dom.nz,j,i);
        fluxLimits(idRW,1,dom.nz,j,i) = fluxLimits(idRW,0,dom.nz,j,i);
        fluxLimits(idTH,1,dom.nz,j,i) = fluxLimits(idTH,0,dom.nz,j,i);
      }
    }

    // Local lax-friedrichs fluxes
    for (int k=0; k<dom.nz+1; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState> s1, s2, f1, f2, upw;
          for (int l=0; l<numState; l++) {
            s1(l) = stateLimits(l,0,k,j,i);
            s2(l) = stateLimits(l,1,k,j,i);
            f1(l) = fluxLimits (l,0,k,j,i);
            f2(l) = fluxLimits (l,1,k,j,i);
          }
          riem.riemannZ(s1, s2, f1, f2, upw);
          for (int l=0; l<numState; l++) {
            flux(l,k,j,i) = upw(l);
          }
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx; i++) {
            tend(l,k,j,i) = - ( flux(l,k+1,j,i) - flux(l,k,j,i) ) / dom.dz;
          }
        }
      }
    }
  }


  inline void compEulerTendSD_S(Array<real> &state, Domain &dom, Array<real> &tend) {
    // Form the tendencies
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(idR ,k,j,i) = 0;
          tend(idRU,k,j,i) = 0;
          tend(idRV,k,j,i) = 0;
          tend(idRW,k,j,i) = -state(idR,hs+k,hs+j,hs+i) * GRAV;
          tend(idTH,k,j,i) = 0;
        }
      }
    }
  }



};

#endif
