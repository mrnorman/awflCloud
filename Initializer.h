
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Hydrostasis.h"
#include "TransformMatrices.h"
#include "Array.h"

class Initializer{

public:

  inline void initialize(State &state, Domain &dom, Parallel &par) {
    SArray<real,ord> gllOrdPoints;
    SArray<real,ord> gllOrdWeights;
    SArray<real,tord> gllTordPoints;
    SArray<real,tord> gllTordWeights;

    //Get GLL points and weights
    TransformMatrices<real> trans;
    trans.get_gll_points(gllOrdPoints);
    trans.get_gll_weights(gllOrdWeights);
    trans.get_gll_points(gllTordPoints);
    trans.get_gll_weights(gllTordWeights);

    // Mock MPI init
    par.nTasks = 1;
    par.myTask = 0;
    par.ix = 0;
    par.iy = 0;
    par.masterTask = 1;
    dom.nx = dom.nx_glob;
    dom.ny = dom.ny_glob;

    // Initialize the grid
    dom.nz = dom.nz_glob;
    dom.dx = dom.xlen / dom.nx_glob;
    dom.dy = dom.ylen / dom.ny_glob;
    dom.dz = dom.zlen / dom.nz_glob;

    // Allocate the fluid state variable
    state.state.setup( numState , dom.nz+2*hs , dom.ny+2*hs , dom.nx+2*hs );

    state.hyDensCells     .setup( dom.nz );
    state.hyDensThetaCells.setup( dom.nz );
    state.hyPressureCells .setup( dom.nz );

    state.hyDensGLL     .setup( dom.nz , tord );
    state.hyDensThetaGLL.setup( dom.nz , tord );
    state.hyPressureGLL .setup( dom.nz , tord );

    Hydrostasis hydro;

    // Initialize the hydrostatic background state for cell averages
    for (int k=0; k<dom.nz; k++) {
      state.hyDensCells     (k) = 0;
      state.hyDensThetaCells(k) = 0;
      state.hyPressureCells (k) = 0;
      // Perform ord-point GLL quadrature for the cell averages
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
        real const t0 = 300._fp;
        real r, t;

        hydro.hydroConstTheta( t0 , zloc , r );
        t = t0;

        state.hyDensCells     (k) += gllOrdWeights(kk) * r;
        state.hyDensThetaCells(k) += gllOrdWeights(kk) * r*t;
        state.hyPressureCells (k) += gllOrdWeights(kk) * mypow( r*t , GAMMA );
      }
    }

    // Initialize the hydrostatic background state for GLL points
    for (int k=0; k<dom.nz; k++) {
      // Perform ord-point GLL quadrature for the cell averages
      for (int kk=0; kk<tord; kk++) {
        real zloc = (k + 0.5_fp)*dom.dz + gllTordPoints(kk)*dom.dz;
        real const t0 = 300._fp;
        real r, t;

        hydro.hydroConstTheta( t0 , zloc , r );
        t = t0;

        state.hyDensGLL     (k,kk) = r;
        state.hyDensThetaGLL(k,kk) = r*t;
        state.hyPressureGLL (k,kk) = mypow( r*t , GAMMA );
      }
    }

    // Initialize the state
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          // Initialize the state to zero
          for (int l=0; l<numState; l++) {
            state.state(l,hs+k,hs+j,hs+i) = 0;
          }
          // Perform ord-point GLL quadrature for the cell averages
          for (int kk=0; kk<ord; kk++) {
            for (int jj=0; jj<ord; jj++) {
              for (int ii=0; ii<ord; ii++) {
                real xloc = (i + 0.5_fp)*dom.dx + gllOrdPoints(ii)*dom.dx;
                real yloc = (j + 0.5_fp)*dom.dy + gllOrdPoints(jj)*dom.dy;
                real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
                real const t0 = 300._fp;
                real r, t;

                hydro.hydroConstTheta( t0 , zloc , r );
                t = t0;
                t += ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 2000, 2000, 2);

                real wt = gllOrdWeights(ii)*gllOrdWeights(jj)*gllOrdWeights(kk);
                state.state(idR ,hs+k,hs+j,hs+i) += wt * r;
                state.state(idRT,hs+k,hs+j,hs+i) += wt * r*t;
              }
            }
          }
        }
      }
    }

    dom.dt = 1.e12_fp;
    // Compute the time step based on the CFL value
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          // Grab state variables
          real r = state.state(idR ,hs+k,hs+j,hs+i)    ;
          real u = state.state(idRU,hs+k,hs+j,hs+i) / r;
          real v = state.state(idRV,hs+k,hs+j,hs+i) / r;
          real w = state.state(idRW,hs+k,hs+j,hs+i) / r;
          real t = state.state(idRT,hs+k,hs+j,hs+i) / r;
          real p = C0 * mypow( r*t , GAMMA );
          real cs = mysqrt( GAMMA * p / r );

          // Compute the max wave
          real maxWave = max( max( myfabs(u) , myfabs(v)) , myfabs(w)) + cs;

          // Compute the time step
          dom.dt = min( dom.dt , dom.cfl * dom.dx / maxWave );
        }
      }
    }

    std::cout << "dx: " << dom.dx << "\n";
    std::cout << "dy: " << dom.dy << "\n";
    std::cout << "dz: " << dom.dz << "\n";
    std::cout << "dt: " << dom.dt << "\n";

  }


  inline real ellipsoid_linear(real const x   , real const y   , real const z ,
                               real const x0  , real const y0  , real const z0,
                               real const xrad, real const yrad, real const zrad, real const amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = mysqrt( xn*xn + yn*yn + zn*zn );
    return amp * max( 1._fp - dist , 0._fp );
  }

};

#endif
