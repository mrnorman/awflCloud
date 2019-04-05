
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Array.h"
#include "Domain.h"
#include "TransformMatrices.h"

class Tendencies {

  Array<real> stateLimits;
  Array<real> fluxLimits;
  Array<real> flux;
  TransformMatrices<real> trans;

public :

  inline void initialize(Domain &dom) {
    fluxLimits .setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    stateLimits.setup(numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    flux       .setup(numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
  }

  inline void compEulerTendSemiX(Array<real> &state, Array<real> &tend, Domain &dom, Parallel &par) {
    SArray<real,ord,ord,ord> s2g_lower_tmp;
    SArray<real,ord,2> s2g_lower;

    // Setup the matrix to transform a stencil of ord cell averages into 2 GLL points
    trans.sten_to_gll_lower( s2g_lower_tmp );
    for (int j=0; j<ord; j++) {
      for (int i=0; i<2; i++) {
        s2g_lower(j,i) = s2g_lower_tmp(1,j,i);
      }
    }

    // Reconstruct to 2 GLL points in the x-direction
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          SArray<real,numState,2> gllState;
          SArray<real,numState,2> gllFlux;

          // Compute GLL points from cell averages
          for (int l=0; l<numState l++) {
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

            fluxGLL(idR ,ii) = r*u;
            fluxGLL(idRU,ii) = r*u*u + p;
            fluxGLL(idRV,ii) = r*u*v;
            fluxGLL(idRW,ii) = r*u*w;
            fluxGLL(idRT,ii) = r*u*t;
          }

          // Store state and flux limits into a globally indexed array
          for (int l=0; l<numState; l++) {
            // Store the left cell edge state and flux estimates
            stateLimits(numState,1,k,j,i  ) = stateGLL(l,0);
            fluxLimits (numState,1,k,j,i  ) = fluxGLL (l,0);

            // Store the Right cell edge state and flux estimates
            stateLimits(numState,0,k,j,i+1) = stateGLL(l,1);
            fluxLimits (numState,0,k,j,i+1) = fluxGLL (l,1);
          }

        }
      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    haloInit      ();
    edgePackN_x   (dom, stateLimits, numState);
    edgePackN_x   (dom, fluxLimits , numState);
    edgeExchange_x(dom, par);
    edgeUnpackN_x (dom, stateLimits, numState);
    edgeUnpackN_x (dom, fluxLimits , numState);

    // Local lax-friedrichs fluxes
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx+1; i++) {
          real r = 0.5_fp * ( stateLimits(idR ,0,k,j,i) + stateLimits(idR ,1,k,j,i) );
          real u = 0.5_fp * ( stateLimits(idRU,0,k,j,i) + stateLimits(idRU,1,k,j,i) ) / r;
          real t = 0.5_fp * ( stateLimits(idRT,0,k,j,i) + stateLimits(idRT,1,k,j,i) ) / r;
          real p = C0 * mypow( r*t , GAMMA );
          real cs = mysqrt( GAMMA * p / r );
          real maxwave = myabs(u) + cs;

          for (int l=0; l<numState; l++) {
            flux(k,j,i) = 0.5_fp * ( fluxLimits(l,1,k,j,i) + fluxLimits(l,0,k,j,i) - maxwave * ( stateLimits(l,1,k,j,i) - stateLimits(l,0,k,j,i) ) )
          }
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int i=0; i<dom.nx+1; i++) {
            tend(k,j,i) = - ( flux(k,j,i+1) - flux(k,j,i) ) / dx;
          }
        }
      }
    }

  }

};

#endif
