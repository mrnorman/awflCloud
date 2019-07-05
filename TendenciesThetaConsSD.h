
#ifndef _TENDENCIES_THETA_CONS_SD_H_
#define _TENDENCIES_THETA_CONS_SD_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
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
  SArray<real,tord> gllWts;
  SArray<real,ord,tord> to_gll;
  SArray<real,ord,ord,ord> wenoRecon;
  SArray<real,hs+2> wenoIdl;
  real wenoSigma;

public :


  inline void initialize(Domain const &dom) {
    TransformMatrices<real> trans;

    fluxLimits  = real5d("fluxLimits",numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    stateLimits = real5d("srcLimits" ,numState,2,dom.nz+1,dom.ny+1,dom.nx+1);
    flux        = real4d("flux"      ,numState  ,dom.nz+1,dom.ny+1,dom.nx+1);
    src         = real3d("src"       ,dom.nz,dom.ny,dom.nx);
    stateGLL    = real5d("stateGLL"  ,numState,dom.nz,dom.ny,dom.nx,tord);

    SArray<real,ord,ord,ord> to_gll_tmp;

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    if (dom.doWeno) {
      trans.coefs_to_gll_lower( to_gll_tmp );
    } else {
      trans.sten_to_gll_lower( to_gll_tmp );
    }
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        to_gll(j,i) = to_gll_tmp(tord-1,j,i);
      }
    }

    trans.weno_sten_to_coefs(wenoRecon);

    trans.get_gll_weights(gllWts);

    wenoSetIdealSigma(wenoIdl,wenoSigma);

  }


  // Transform ord stencil cell averages into tord GLL point values
  inline _HOSTDEV void reconStencil(SArray<real,ord> const &stencil, SArray<real,tord> &gll, int const doWeno,
                                    SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll,
                                    SArray<real,hs+2> const &wenoIdl, real wenoSigma) {
    SArray<real,ord> coefs;
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


  inline void compEulerTend_X(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {


    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    reconSD_X(state, dom.hyDensCells, dom.hyDensThetaCells, dom, wenoRecon, to_gll, stateLimits, fluxLimits, wenoIdl, wenoSigma);

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_x   (dom, stateLimits, numState);
    exch.edgePackN_x   (dom, fluxLimits , numState);
    exch.edgeExchange_x(dom, par);
    exch.edgeUnpackN_x (dom, stateLimits, numState);
    exch.edgeUnpackN_x (dom, fluxLimits , numState);

    // Riemann solver
    computeFlux_X(stateLimits, fluxLimits, flux, dom);

    // Form the tendencies
    computeTend_X(flux, tend, dom);
  }


  inline void compEulerTend_Y(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {

    //Exchange halos in the y-direction
    exch.haloInit      ();
    exch.haloPackN_y   (dom, state, numState);
    exch.haloExchange_y(dom, par);
    exch.haloUnpackN_y (dom, state, numState);

    // Reconstruct to tord GLL points in the y-direction
    reconSD_Y(state, dom.hyDensCells, dom.hyDensThetaCells, dom, wenoRecon, to_gll, stateLimits, fluxLimits, wenoIdl, wenoSigma);

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_y   (dom, stateLimits, numState);
    exch.edgePackN_y   (dom, fluxLimits , numState);
    exch.edgeExchange_y(dom, par);
    exch.edgeUnpackN_y (dom, stateLimits, numState);
    exch.edgeUnpackN_y (dom, fluxLimits , numState);

    // Riemann solver
    computeFlux_Y(stateLimits, fluxLimits, flux, dom);

    // Form the tendencies
    computeTend_Y(flux, tend, dom);
  }


  inline void compEulerTend_Z(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {

    // Boundaries for the fluid state in the z-direction
    stateBoundariesZ(state, dom);

    // Reconstruct to tord GLL points in the x-direction
    reconSD_Z(state, dom.hyDensGLL, dom.hyDensThetaGLL, dom, wenoRecon, to_gll, stateLimits, fluxLimits, wenoIdl, wenoSigma);

    // Apply boundary conditions to fluxes and state values
    edgeBoundariesZ(stateLimits, fluxLimits, dom);

    // Riemann solver
    computeFlux_Z(stateLimits, fluxLimits, flux, dom);

    // Form the tendencies
    computeTend_Z(flux, tend, src, dom);
  }


  inline void compEulerTend_S(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    // Form the tendencies
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      tend(idR ,k,j,i) = 0;
      tend(idRU,k,j,i) = 0;
      tend(idRV,k,j,i) = 0;
      tend(idRW,k,j,i) = -state(idR,hs+k,hs+j,hs+i) * GRAV;
      tend(idRT,k,j,i) = 0;
    });
  }


  inline void reconSD_X(real4d &state, real1d const &hyDensCells, real1d const &hyDensThetaCells,
                        Domain const &dom, SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll, 
                        real5d &stateLimits, real5d &fluxLimits, SArray<real,hs+2> const &wenoIdl, real &wenoSigma) {
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      SArray<real,numState,tord> gllState;  // GLL state values
      SArray<real,numState,tord> gllFlux;   // GLL flux values

      // Compute tord GLL points of the state vector
      for (int l=0; l<numState; l++) {
        SArray<real,ord> stencil;
        SArray<real,tord> gllPts;
        for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,hs+j,i+ii); }
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
      }
      for (int ii=0; ii<tord; ii++) {
        gllState(idR ,ii) += dom.hyDensCells     (hs+k);
        gllState(idRT,ii) += dom.hyDensThetaCells(hs+k);
      }

      // Compute fluxes and at the GLL points
      for (int ii=0; ii<tord; ii++) {
        real r = gllState(idR ,ii);
        real u = gllState(idRU,ii) / r;
        real v = gllState(idRV,ii) / r;
        real w = gllState(idRW,ii) / r;
        real t = gllState(idRT,ii) / r;
        real p = C0 * pow( r*t , GAMMA );

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
        stateLimits(l,0,k,j,i+1) = gllState(l,tord-1);
        fluxLimits (l,0,k,j,i+1) = gllFlux (l,tord-1);
      }

    });
  }


  inline void reconSD_Y(real4d &state, real1d const &hyDensCells, real1d const &hyDensThetaCells,
                        Domain const &dom, SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll, 
                        real5d &stateLimits, real5d &fluxLimits, SArray<real,hs+2> const &wenoIdl, real &wenoSigma) {
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      SArray<real,numState,tord> gllState;
      SArray<real,numState,tord> gllFlux;

      // Compute GLL points from cell averages
      for (int l=0; l<numState; l++) {
        SArray<real,ord> stencil;
        SArray<real,tord> gllPts;
        for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+k,j+ii,hs+i); }
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
      }
      for (int ii=0; ii<tord; ii++) {
        gllState(idR ,ii) += hyDensCells     (hs+k);
        gllState(idRT,ii) += hyDensThetaCells(hs+k);
      }

      // Compute fluxes and at the GLL points
      for (int ii=0; ii<tord; ii++) {
        real r = gllState(idR ,ii);
        real u = gllState(idRU,ii) / r;
        real v = gllState(idRV,ii) / r;
        real w = gllState(idRW,ii) / r;
        real t = gllState(idRT,ii) / r;
        real p = C0 * pow( r*t , GAMMA );

        gllFlux(idR ,ii) = r*v;
        gllFlux(idRU,ii) = r*v*u;
        gllFlux(idRV,ii) = r*v*v + p;
        gllFlux(idRW,ii) = r*v*w;
        gllFlux(idRT,ii) = r*v*t;
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

    });
  }


  inline void reconSD_Z(real4d &state, real2d const &hyDensGLL, real2d const &hyDensThetaGLL,
                        Domain const &dom, SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll, 
                        real5d &stateLimits, real5d &fluxLimits, SArray<real,hs+2> const &wenoIdl, real &wenoSigma) {
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      SArray<real,numState,tord> gllState;
      SArray<real,numState,tord> gllFlux;

      // Compute GLL points from cell averages
      for (int l=0; l<numState; l++) {
        SArray<real,ord> stencil;
        SArray<real,tord> gllPts;
        for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,k+ii,hs+j,hs+i); }
        reconStencil(stencil, gllPts, dom.doWeno, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
      }
      for (int ii=0; ii<tord; ii++) {
        gllState(idR ,ii) += hyDensGLL     (k,ii);
        gllState(idRT,ii) += hyDensThetaGLL(k,ii);
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
        real t = gllState(idRT,ii) / r;
        real p = C0 * pow( r*t , GAMMA );

        gllFlux(idR ,ii) = r*w;
        gllFlux(idRU,ii) = r*w*u;
        gllFlux(idRV,ii) = r*w*v;
        gllFlux(idRW,ii) = r*w*w + p - C0*pow(hyDensThetaGLL(k,ii),GAMMA);
        gllFlux(idRT,ii) = r*w*t;
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

    });
  }


  inline void computeFlux_X(real5d const &stateLimits, real5d const &fluxLimits, real4d &flux, Domain const &dom ) {
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx+1; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*(dom.nx+1) , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx+1,k,j,i);
      SArray<real,numState> s1, s2, f1, f2, upw;
      for (int l=0; l<numState; l++) {
        s1(l) = stateLimits(l,0,k,j,i);
        s2(l) = stateLimits(l,1,k,j,i);
        f1(l) = fluxLimits (l,0,k,j,i);
        f2(l) = fluxLimits (l,1,k,j,i);
      }
      riemannX(s1, s2, f1, f2, upw);
      for (int l=0; l<numState; l++) {
        flux(l,k,j,i) = upw(l);
        // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dx/dom.dt * (s2(l) - s1(l)) );
      }
    });
  }


  inline void computeTend_X(real4d const &flux, real4d &tend, Domain const &dom) {
    // for (int l=0; l<numState; l++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    //       for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( numState*dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int l, k, j, i;
      unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
      tend(l,k,j,i) = - ( flux(l,k,j,i+1) - flux(l,k,j,i) ) / dom.dx;
    });
  }


  inline void computeFlux_Y(real5d const &stateLimits, real5d const &fluxLimits, real4d &flux, Domain const &dom ) {
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny+1; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*(dom.ny+1)*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny+1,dom.nx,k,j,i);
      SArray<real,numState> s1, s2, f1, f2, upw;
      for (int l=0; l<numState; l++) {
        s1(l) = stateLimits(l,0,k,j,i);
        s2(l) = stateLimits(l,1,k,j,i);
        f1(l) = fluxLimits (l,0,k,j,i);
        f2(l) = fluxLimits (l,1,k,j,i);
      }
      riemannY(s1, s2, f1, f2, upw);
      for (int l=0; l<numState; l++) {
        flux(l,k,j,i) = upw(l);
        // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dy/dom.dt * (s2(l) - s1(l)) );
      }
    });
  }


  inline void computeTend_Y(real4d const &flux, real4d &tend, Domain const &dom) {
    // for (int l=0; l<numState; l++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    //       for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( numState*dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int l, k, j, i;
      unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
      tend(l,k,j,i) = - ( flux(l,k,j+1,i) - flux(l,k,j,i) ) / dom.dy;
    });
  }


  inline void computeFlux_Z(real5d const &stateLimits, real5d const &fluxLimits, real4d &flux, Domain const &dom ) {
    // for (int k=0; k<dom.nz+1; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( (dom.nz+1)*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz+1,dom.ny,dom.nx,k,j,i);
      SArray<real,numState> s1, s2, f1, f2, upw;
      for (int l=0; l<numState; l++) {
        s1(l) = stateLimits(l,0,k,j,i);
        s2(l) = stateLimits(l,1,k,j,i);
        f1(l) = fluxLimits (l,0,k,j,i);
        f2(l) = fluxLimits (l,1,k,j,i);
      }
      riemannZ(s1, s2, f1, f2, upw);
      for (int l=0; l<numState; l++) {
        flux(l,k,j,i) = upw(l);
        // flux(l,k,j,i) = 0.5_fp * ( f2(l) + f1(l) - dom.cfl*dom.dz/dom.dt * (s2(l) - s1(l)) );
      }
    });
  }


  inline void computeTend_Z(real4d const &flux, real4d &tend, real3d &src, Domain const &dom) {
    // for (int l=0; l<numState; l++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    //       for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( numState*dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int l, k, j, i;
      unpackIndices(iGlob,numState,dom.nz,dom.ny,dom.nx,l,k,j,i);
      tend(l,k,j,i) = - ( flux(l,k+1,j,i) - flux(l,k,j,i) ) / dom.dz;
      if (l==idRW) {
        tend(l,k,j,i) += src(k,j,i);
      }
    });
  }


  inline void stateBoundariesZ(real4d &state, Domain const &dom) {
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    //     for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( dom.ny*dom.nx*hs , KOKKOS_LAMBDA (int const iGlob) {
      int j, i, ii;
      unpackIndices(iGlob,dom.ny,dom.nx,hs,j,i,ii);
      state(idR ,ii,hs+j,hs+i) = state(idR ,hs,hs+j,hs+i);
      state(idRU,ii,hs+j,hs+i) = state(idRU,hs,hs+j,hs+i);
      state(idRV,ii,hs+j,hs+i) = state(idRV,hs,hs+j,hs+i);
      state(idRW,ii,hs+j,hs+i) = 0;
      state(idRT,ii,hs+j,hs+i) = state(idRT,hs,hs+j,hs+i);

      state(idR ,dom.nz+hs+ii,hs+j,hs+i) = state(idR ,dom.nz+hs-1,hs+j,hs+i);
      state(idRU,dom.nz+hs+ii,hs+j,hs+i) = state(idRU,dom.nz+hs-1,hs+j,hs+i);
      state(idRV,dom.nz+hs+ii,hs+j,hs+i) = state(idRV,dom.nz+hs-1,hs+j,hs+i);
      state(idRW,dom.nz+hs+ii,hs+j,hs+i) = 0;
      state(idRT,dom.nz+hs+ii,hs+j,hs+i) = state(idRT,dom.nz+hs-1,hs+j,hs+i);
    });
  }


  inline void edgeBoundariesZ(real5d &stateLimits, real5d &fluxLimits, Domain const &dom) {
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int j, i;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      stateLimits(idR ,0,0     ,j,i) = stateLimits(idR ,1,0     ,j,i);
      stateLimits(idRU,0,0     ,j,i) = stateLimits(idRU,1,0     ,j,i);
      stateLimits(idRV,0,0     ,j,i) = stateLimits(idRV,1,0     ,j,i);
      stateLimits(idRW,0,0     ,j,i) = 0;
      stateLimits(idRW,1,0     ,j,i) = 0;
      stateLimits(idRT,0,0     ,j,i) = stateLimits(idRT,1,0     ,j,i);

      stateLimits(idR ,1,dom.nz,j,i) = stateLimits(idR ,0,dom.nz,j,i);
      stateLimits(idRU,1,dom.nz,j,i) = stateLimits(idRU,0,dom.nz,j,i);
      stateLimits(idRV,1,dom.nz,j,i) = stateLimits(idRV,0,dom.nz,j,i);
      stateLimits(idRW,0,dom.nz,j,i) = 0;
      stateLimits(idRW,1,dom.nz,j,i) = 0;
      stateLimits(idRT,1,dom.nz,j,i) = stateLimits(idRT,0,dom.nz,j,i);

      fluxLimits(idR ,0,0     ,j,i) = fluxLimits(idR ,1,0     ,j,i);
      fluxLimits(idRU,0,0     ,j,i) = fluxLimits(idRU,1,0     ,j,i);
      fluxLimits(idRV,0,0     ,j,i) = fluxLimits(idRV,1,0     ,j,i);
      fluxLimits(idRW,0,0     ,j,i) = fluxLimits(idRW,1,0     ,j,i);
      fluxLimits(idRT,0,0     ,j,i) = fluxLimits(idRT,1,0     ,j,i);

      fluxLimits(idR ,1,dom.nz,j,i) = fluxLimits(idR ,0,dom.nz,j,i);
      fluxLimits(idRU,1,dom.nz,j,i) = fluxLimits(idRU,0,dom.nz,j,i);
      fluxLimits(idRV,1,dom.nz,j,i) = fluxLimits(idRV,0,dom.nz,j,i);
      fluxLimits(idRW,1,dom.nz,j,i) = fluxLimits(idRW,0,dom.nz,j,i);
      fluxLimits(idRT,1,dom.nz,j,i) = fluxLimits(idRT,0,dom.nz,j,i);
    });
  }


  inline void compStrakaTend(real4d &state, Domain const &dom, Exchange &exch, Parallel const &par, real4d &tend) {
    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    //Exchange halos in the y-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);
    
    // Boundaries for the fluid state in the z-direction
    stateBoundariesZ(state, dom);

    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      real r = ( state(idR,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k) );
      SArray<real,numState,3> sten;

      tend(idR,k,j,i) = 0.;

      real dx2 = dom.dx*dom.dx;
      for (int l=1; l<numState; l++) {
        for (int ii=-1; ii<=1; ii++) {
          if (l == idRT) {
            sten(l,ii+1) = ( state(l,hs+k,hs+j,hs+i+ii) + dom.hyDensThetaCells(hs+k) ) / ( state(idR,hs+k,hs+j,hs+i+ii) + dom.hyDensCells(hs+k) );
          } else {
            sten(l,ii+1) = state(l,hs+k,hs+j,hs+i+ii) / ( state(idR,hs+k,hs+j,hs+i+ii) + dom.hyDensCells(hs+k) );
          }
        }
        tend(l,k,j,i)  = 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dx2;
      }

      real dy2 = dom.dy*dom.dy;
      for (int l=1; l<numState; l++) {
        for (int ii=-1; ii<=1; ii++) {
          if (l == idRT) {
            sten(l,ii+1) = ( state(l,hs+k,hs+j+ii,hs+i) + dom.hyDensThetaCells(hs+k) ) / ( state(idR,hs+k,hs+j+ii,hs+i) + dom.hyDensCells(hs+k) );
          } else {
            sten(l,ii+1) = state(l,hs+k,hs+j+ii,hs+i) / ( state(idR,hs+k,hs+j+ii,hs+i) + dom.hyDensCells(hs+k) );
          }
        }
        tend(l,k,j,i) += 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dy2;
      }

      real dz2 = dom.dz*dom.dz;
      for (int l=1; l<numState; l++) {
        for (int ii=-1; ii<=1; ii++) {
          if (l == idRT) {
            sten(l,ii+1) = ( state(l,hs+k+ii,hs+j,hs+i) + dom.hyDensThetaCells(hs+k+ii) ) / ( state(idR,hs+k+ii,hs+j,hs+i) + dom.hyDensCells(hs+k+ii) );
          } else {
            sten(l,ii+1) = state(l,hs+k+ii,hs+j,hs+i) / ( state(idR,hs+k+ii,hs+j,hs+i) + dom.hyDensCells(hs+k+ii) );
          }
        }
        tend(l,k,j,i) += 75 * r * ( sten(l,2) - 2*sten(l,1) + sten(l,0) ) / dz2;
      }
    });
  }


};

#endif
