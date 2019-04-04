
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

    state.hyDensCells     .setup( dom.nz+2*hs );
    state.hyDensThetaCells.setup( dom.nz+2*hs );
    state.hyPressureCells .setup( dom.nz+2*hs );

    state.hyDensGLL     .setup( dom.nz , tord );
    state.hyDensThetaGLL.setup( dom.nz , tord );
    state.hyPressureGLL .setup( dom.nz , tord );

    Hydrostasis hydro;

    // Initialize the hydrostatic background state for cell averages
    for (int k=0; k<dom.nz; k++) {
      state.hyDensCells     (hs+k) = 0;
      state.hyDensThetaCells(hs+k) = 0;
      state.hyPressureCells (hs+k) = 0;
      // Perform ord-point GLL quadrature for the cell averages
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
        real const t0 = 300._fp;
        real r, t;

        hydro.hydroConstTheta( t0 , zloc , r );
        t = t0;

        state.hyDensCells     (hs+k) += gllOrdWeights(kk) * r;
        state.hyDensThetaCells(hs+k) += gllOrdWeights(kk) * r*t;
        state.hyPressureCells (hs+k) += gllOrdWeights(kk) * mypow( r*t , GAMMA );
      }
    }
    //TODO: Apply boundary conditions to background state for cell averages

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

    // Compute the time step based on the CFL value

  }

};

#endif
