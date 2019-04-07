
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Hydrostasis.h"
#include "Exchange.h"
#include "TransformMatrices.h"
#include "TimeIntegrator.h"
#include "Array.h"
#include "mpi.h"

class Initializer{

public:


  void initialize_mpi( int *argc , char ***argv , Parallel &par ) {
    int ierr = MPI_Init( argc , argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    //Determine if I'm the master process
    if (par.myrank == 0) {
      par.masterproc = 1;
    } else {
      par.masterproc = 0;
    }
  }

  void initialize(State &state, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint) {
    int ierr;
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
    if (par.nranks != par.nproc_x*par.nproc_y) {
      std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
      std::cerr << par.nproc_x << " " << par.nproc_y << " " << par.nranks << "\n";
      exit(-1);
    }

    //Get my x and y process grid ID
    par.px = par.myrank % par.nproc_x;
    par.py = par.myrank / par.nproc_x;

    //Get my beginning and ending global indices
    double nper;
    nper = ((double) dom.nx_glob)/par.nproc_x;
    par.i_beg = (long) round( nper* par.px    );
    par.i_end = (long) round( nper*(par.px+1) )-1;
    nper = ((double) dom.ny_glob)/par.nproc_y;
    par.j_beg = (long) round( nper* par.py    );
    par.j_end = (long) round( nper*(par.py+1) )-1;
    //Determine my number of grid cells
    dom.nx = par.i_end - par.i_beg + 1;
    dom.ny = par.j_end - par.j_beg + 1;
    par.neigh.setup(3,3);
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int pxloc = par.px+i-1;
        if (pxloc < 0            ) pxloc = pxloc + par.nproc_x;
        if (pxloc > par.nproc_x-1) pxloc = pxloc - par.nproc_x;
        int pyloc = par.py+j-1;
        if (pyloc < 0            ) pyloc = pyloc + par.nproc_y;
        if (pyloc > par.nproc_y-1) pyloc = pyloc - par.nproc_y;
        par.neigh(j,i) = pyloc * par.nproc_x + pxloc;
      }
    }

    // Debug output for the parallel decomposition
    if (0) {
      for (int rr=0; rr < par.nranks; rr++) {
        if (rr == par.myrank) {
          std::cout << "Hello! My Rank is what, my rank is who, my rank is: " << par.myrank << "\n";
          std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
          std::cout << "I have: " << dom.nx << " x " << dom.ny << " columns." << "\n";
          std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
          std::cout << "My neighbor matrix is:\n";
          for (int j = 2; j >= 0; j--) {
            for (int i = 0; i < 3; i++) {
              std::cout << std::setw(6) << par.neigh(j,i) << " ";
            }
            printf("\n");
          }
          printf("\n");
        }
        ierr = MPI_Barrier(MPI_COMM_WORLD);
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }

    // Initialize the grid
    dom.etime = 0;
    dom.nz = dom.nz_glob;
    dom.dx = dom.xlen / dom.nx_glob;
    dom.dy = dom.ylen / dom.ny_glob;
    dom.dz = dom.zlen / dom.nz_glob;

    tint.initialize(dom);

    // Allocate the MPI exchange buffers
    exch.allocate(dom);

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

    // Enforce vertical boundaries
    for (int ii=0; ii<hs; ii++) {
      state.hyDensCells     (ii) = state.hyDensCells     (hs);
      state.hyDensThetaCells(ii) = state.hyDensThetaCells(hs);
      state.hyPressureCells (ii) = state.hyPressureCells (hs);
      state.hyDensCells     (dom.nz+hs+ii) = state.hyDensCells     (dom.nz+hs-1);
      state.hyDensThetaCells(dom.nz+hs+ii) = state.hyDensThetaCells(dom.nz+hs-1);
      state.hyPressureCells (dom.nz+hs+ii) = state.hyPressureCells (dom.nz+hs-1);
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
                real xloc = (par.i_beg + i + 0.5_fp)*dom.dx + gllOrdPoints(ii)*dom.dx;
                real yloc = (par.j_beg + j + 0.5_fp)*dom.dy + gllOrdPoints(jj)*dom.dy;
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

    real dtloc = dom.dt;
    ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_REAL , MPI_MIN, MPI_COMM_WORLD);


    if (par.masterproc) {
      std::cout << "dx: " << dom.dx << "\n";
      std::cout << "dy: " << dom.dy << "\n";
      std::cout << "dz: " << dom.dz << "\n";
      std::cout << "dt: " << dom.dt << "\n";
    }

  }


  real ellipsoid_linear(real const x   , real const y   , real const z ,
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
