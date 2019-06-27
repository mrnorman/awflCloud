
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Hydrostasis.h"
#include "Exchange.h"
#include "TransformMatrices.h"
#include "TimeIntegrator.h"
#include "Indexing.h"
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

    if (dom.ny_glob == 1) { dom.run2d = 1; } else { dom.run2d = 0; }

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
    state.state = real4d( "state" , numState , dom.nz+2*hs , dom.ny+2*hs , dom.nx+2*hs );

    state.hyDensCells      = real1d( "hyCellsR"  , dom.nz+2*hs );
    state.hyDensThetaCells = real1d( "hyCellsRT" , dom.nz+2*hs );
    state.hyThetaCells     = real1d( "hyCellsT"  , dom.nz+2*hs );
    state.hyPressureCells  = real1d( "hyCellsp"  , dom.nz+2*hs );

    state.hyDensGLL      = real2d( "hyGLLR"  , dom.nz , tord );
    state.hyDensThetaGLL = real2d( "hyGLLRT" , dom.nz , tord );
    state.hyThetaGLL     = real2d( "hyGLLT"  , dom.nz , tord );
    state.hyPressureGLL  = real2d( "hyGLLp"  , dom.nz , tord );

    // Initialize the hydrostatic background state for cell averages
    // for (int k=0; k<dom.nz; k++) {
    Kokkos::parallel_for( "initHydroCells" , dom.nz , KOKKOS_LAMBDA (int const k) {
      state.hyDensCells     (hs+k) = 0;
      state.hyDensThetaCells(hs+k) = 0;
      state.hyPressureCells (hs+k) = 0;
      // Perform ord-point GLL quadrature for the cell averages
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
        real r0, t0;

        if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
          t0 = 300._fp;
          hydro::hydroConstTheta( t0 , zloc , r0 );
        }

        state.hyDensCells     (hs+k) += gllOrdWeights(kk) * r0;
        state.hyDensThetaCells(hs+k) += gllOrdWeights(kk) * r0*t0;
        state.hyThetaCells    (hs+k) += gllOrdWeights(kk) * t0;
        state.hyPressureCells (hs+k) += gllOrdWeights(kk) * pow( r0*t0 , GAMMA );
      }
    });

    // Enforce vertical boundaries
    // for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( "boundariesHydroCells" , hs , KOKKOS_LAMBDA (int const ii) {
      state.hyDensCells     (ii) = state.hyDensCells     (hs);
      state.hyDensThetaCells(ii) = state.hyDensThetaCells(hs);
      state.hyThetaCells    (ii) = state.hyThetaCells    (hs);
      state.hyPressureCells (ii) = state.hyPressureCells (hs);
      state.hyDensCells     (dom.nz+hs+ii) = state.hyDensCells     (dom.nz+hs-1);
      state.hyDensThetaCells(dom.nz+hs+ii) = state.hyDensThetaCells(dom.nz+hs-1);
      state.hyThetaCells    (dom.nz+hs+ii) = state.hyThetaCells    (dom.nz+hs-1);
      state.hyPressureCells (dom.nz+hs+ii) = state.hyPressureCells (dom.nz+hs-1);
    });

    // Initialize the hydrostatic background state for GLL points
    // Perform ord-point GLL quadrature for the cell averages
    // for (int k=0; k<dom.nz; k++) {
    //   for (int kk=0; kk<tord; kk++) {
    Kokkos::parallel_for( "initHydroGLL" , dom.nz*tord , KOKKOS_LAMBDA (int const iGlob) {
      int k, kk;
      unpackIndices(iGlob,dom.nz,tord,k,kk);
      real zloc = (k + 0.5_fp)*dom.dz + gllTordPoints(kk)*dom.dz;
      real r0, t0;

      if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
        t0 = 300._fp;
        hydro::hydroConstTheta( t0 , zloc , r0 );
      }

      state.hyDensGLL     (k,kk) = r0;
      state.hyDensThetaGLL(k,kk) = r0*t0;
      state.hyThetaGLL    (k,kk) = t0;
      state.hyPressureGLL (k,kk) = pow( r0*t0 , GAMMA );
    });

    // Initialize the state
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( "InitFluidState" , dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
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
            real r0, t0, r, t;

            if (dom.run2d) yloc = dom.ylen/2;

            if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
              t0 = 300._fp;
              hydro::hydroConstTheta( t0 , zloc , r0 );
            }

            if        (dom.dataInit == DATA_INIT_COLLISION) {
              t  = ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 2000, 2000,  20);
              t += ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 8000, 2000, 2000, 2000, -20);
            } else if (dom.dataInit == DATA_INIT_THERMAL  ) {
              t  = ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 2000, 2000, 2  );
            } else if (dom.dataInit == DATA_INIT_STRAKA   ) {
              t  = ellipsoid_cosine(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 3000, 4000, 4000, 2000, -15 , 1 );
            }

            // Set the initial density such that pressure is constant (to get rid of distracting acoustic waves)
            r = r0*(t0/(t0+t)-1._fp);

            real wt = gllOrdWeights(ii)*gllOrdWeights(jj)*gllOrdWeights(kk);

            if        (dom.eqnSet == EQN_THETA_CONS) {
              state.state(idR ,hs+k,hs+j,hs+i) += wt * r  ;
              state.state(idRT,hs+k,hs+j,hs+i) += wt * ( (r0+r)*(t0+t) - r0*t0 );
            } else if (dom.eqnSet == EQN_THETA_PRIM) {
              state.state(idR ,hs+k,hs+j,hs+i) += wt * r;
              state.state(idRT,hs+k,hs+j,hs+i) += wt * t;
            }
          }
        }
      }
    });

    real3d dt3d("dt3d",dom.nz,dom.ny,dom.nx);
    // Compute the time step based on the CFL value
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( "Compute_dt3d" , dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      // Grab state variables
      real r = state.state(idR ,hs+k,hs+j,hs+i) + state.hyDensCells(hs+k);
      real u = state.state(idRU,hs+k,hs+j,hs+i) / r;
      real v = state.state(idRV,hs+k,hs+j,hs+i) / r;
      real w = state.state(idRW,hs+k,hs+j,hs+i) / r;
      real t = ( state.state(idRT,hs+k,hs+j,hs+i) + state.hyDensThetaCells(hs+k) ) / r;
      real p = C0 * pow( r*t , GAMMA );
      real cs = sqrt( GAMMA * p / r );

      // Compute the max wave
      real maxWave = max( max( fabs(u) , fabs(v)) , fabs(w)) + cs;

      // Compute the time step
      real mindx = min( min(dom.dx,dom.dy) , dom.dz);
      dt3d(k,j,i) = dom.cfl * mindx / maxWave;
    });

    dom.dt = 1.e12_fp;
    Kokkos::parallel_reduce( "ReduceTimeStep" , dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob, real &dt) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      dt = min(dt,dt3d(k,j,i));
    } , Kokkos::Min<real>(dom.dt) );

    Kokkos::fence();

    real dtloc = dom.dt;
    ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_REAL , MPI_MIN, MPI_COMM_WORLD);

    if (par.masterproc) {
      std::cout << "dx: " << dom.dx << "\n";
      std::cout << "dy: " << dom.dy << "\n";
      std::cout << "dz: " << dom.dz << "\n";
      std::cout << "dt: " << dom.dt << "\n";
    }

  }


  inline _HOSTDEV real ellipsoid_linear(real const x   , real const y   , real const z ,
                                        real const x0  , real const y0  , real const z0,
                                        real const xrad, real const yrad, real const zrad, real const amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    return amp * max( 1._fp - dist , 0._fp );
  }

  inline _HOSTDEV real ellipsoid_cosine(real const x   , real const y   , real const z ,
                                        real const x0  , real const y0  , real const z0,
                                        real const xrad, real const yrad, real const zrad, real const amp, real const pwr) {
    real val = 0;
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    if (dist <= 1._fp) {
      val = amp * pow( (cos(PI*dist)+1)/2 , pwr );
    }
    return val;
  }

};

#endif
