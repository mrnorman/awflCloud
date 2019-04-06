
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Array.h"
#include "Domain.h"
#include "Exchange.h"
#include "TransformMatrices.h"

class Tendencies {

  Array<real> stateLimits;
  Array<real> fluxLimits;
  Array<real> flux;
  TransformMatrices<real> trans;

public :

  Array<real> tend;

  inline void initialize(Domain &dom) {
    fluxLimits .setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    stateLimits.setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    flux       .setup(numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
    tend       .setup(numState  ,dom.nz  ,dom.ny  ,dom.nx  );
  }

  inline void compEulerTendSD_X(Array<real> &state, Domain &dom, Exchange &exch, Parallel &par) {
    SArray<real,ord,ord,ord> s2g_lower_tmp;
    SArray<real,ord,2> s2g_lower;

    // Setup the matrix to transform a stencil of ord cell averages into 2 GLL points
    trans.sten_to_gll_lower( s2g_lower_tmp );
    for (int j=0; j<ord; j++) {
      for (int i=0; i<2; i++) {
        s2g_lower(j,i) = s2g_lower_tmp(1,j,i);
      }
    }

    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    // Reconstruct to 2 GLL points in the x-direction
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState,2> gllState;
          SArray<real,numState,2> gllFlux;

          // Compute GLL points from cell averages
          for (int l=0; l<numState; l++) {
            for (int ii=0; ii<2; ii++) {
              gllState(l,ii) = 0.;
              for (int s=0; s<ord; s++) {
                gllState(l,ii) += s2g_lower(s,ii) * state(l,hs+k,hs+j,i+s);
              }
            }
          }

          // Compute fluxes and at the GLL points
          for (int ii=0; ii<2; ii++) {
            real r = gllState(idR ,ii);
            real u = gllState(idRU,ii) / r;
            real v = gllState(idRV,ii) / r;
            real w = gllState(idRW,ii) / r;
            real t = gllState(idRT,ii) / r;
            real p = C0 * mypow( r*t , GAMMA );

            gllFlux(idR ,ii) = r*u;
            gllFlux(idRU,ii) = r*u*u + p;
            gllFlux(idRV,ii) = r*u*v;
            gllFlux(idRW,ii) = r*u*w;
            gllFlux(idRT,ii) = r*u*t;
          }

          // Store state and flux limits into a globally indexed array
          for (int l=0; l<numState; l++) {
            // Store the left cell edge state and flux estimates
            stateLimits(l,1,k,j,i  ) = gllState(l,0);
            fluxLimits (l,1,k,j,i  ) = gllFlux (l,0);

            // Store the Right cell edge state and flux estimates
            stateLimits(l,0,k,j,i+1) = gllState(l,1);
            fluxLimits (l,0,k,j,i+1) = gllFlux (l,1);
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
          real r = 0.5_fp * ( stateLimits(idR ,0,k,j,i) + stateLimits(idR ,1,k,j,i) );
          real u = 0.5_fp * ( stateLimits(idRU,0,k,j,i) + stateLimits(idRU,1,k,j,i) ) / r;
          real t = 0.5_fp * ( stateLimits(idRT,0,k,j,i) + stateLimits(idRT,1,k,j,i) ) / r;
          real p = C0 * mypow( r*t , GAMMA );
          real cs = mysqrt( GAMMA * p / r );
          real maxwave = myfabs(u) + cs;

          for (int l=0; l<numState; l++) {
            flux(l,k,j,i) = 0.5_fp * ( fluxLimits(l,1,k,j,i) + fluxLimits(l,0,k,j,i) - maxwave * ( stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i) ) );
            // flux(l,k,j,i) = 0.5_fp * ( fluxLimits(l,1,k,j,i) + fluxLimits(l,0,k,j,i)                                                                 );
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

};

#endif
