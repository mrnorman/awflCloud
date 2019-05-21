
#ifndef _FILEIO_H_
#define _FILEIO_H_

#include "const.h"
#include "State.h"
#include "pnetcdf.h"
#include "Indexing.h"
#include "mpi.h"

class FileIO {

protected:

  real outTimer;
  int ncid, numOut;
  int tDim, xDim, yDim, zDim;
  int tVar, xVar, yVar, zVar, rVar, uVar, vVar, wVar, thVar, hyrVar, hyrtVar, pVar;

public:

  void outputInit(State &state, Domain const &dom, Parallel const &par) {
    int dimids[4];
    MPI_Offset st[1], ct[1];
    real1d xCoord("xCoord",dom.nx);
    real1d yCoord("yCoord",dom.ny);
    real1d zCoord("zCoord",dom.nz);

    numOut = 0;

    outTimer = 0.;

    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER | NC_64BIT_DATA , MPI_INFO_NULL , &ncid ) , __LINE__ );
    // ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "z" , (MPI_Offset) dom.nz_glob  , &zDim ) , __LINE__ );

    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_FLOAT , 1 , dimids , &tVar ) , __LINE__ );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_FLOAT , 1 , dimids , &xVar ) , __LINE__ );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_FLOAT , 1 , dimids , &yVar ) , __LINE__ );
    dimids[0] = zDim;
    ncwrap( ncmpi_def_var( ncid , "z"      , NC_FLOAT , 1 , dimids , &zVar ) , __LINE__ );
    dimids[0] = tDim; dimids[1] = zDim; dimids[2] = yDim; dimids[3] = xDim;
    ncwrap( ncmpi_def_var( ncid , "density" , NC_FLOAT , 4 , dimids , &rVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "u"       , NC_FLOAT , 4 , dimids , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "v"       , NC_FLOAT , 4 , dimids , &vVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "w"       , NC_FLOAT , 4 , dimids , &wVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "theta"   , NC_FLOAT , 4 , dimids , &thVar ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "pressure", NC_FLOAT , 4 , dimids , &pVar  ) , __LINE__ );
    dimids[0] = zDim;
    ncwrap( ncmpi_def_var( ncid , "hyDens"      , NC_FLOAT , 1 , dimids , &hyrVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "hyDensTheta" , NC_FLOAT , 1 , dimids , &hyrtVar ) , __LINE__ );

    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Compute x, y, and z coordinates
    // for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nx , KOKKOS_LAMBDA (int const i) {
      xCoord(i) = ( par.i_beg + i + 0.5_fp ) * dom.dx;
    });
    // for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( dom.ny , KOKKOS_LAMBDA (int const j) {
      yCoord(j) = ( par.j_beg + j + 0.5_fp ) * dom.dy;
    });
    // for (int k=0; k<dom.nz; k++) {
    Kokkos::parallel_for( dom.nz , KOKKOS_LAMBDA (int const k) {
      zCoord(k) = (             k + 0.5_fp ) * dom.dz;
    });
    Kokkos::fence();

    // Write out x, y, and z coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap( ncmpi_put_vara_float_all( ncid , xVar , st , ct , xCoord.data() ) , __LINE__ );
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap( ncmpi_put_vara_float_all( ncid , yVar , st , ct , yCoord.data() ) , __LINE__ );

    // Write out the hydrostatic background states and z coordinates
    st[0] = 0;
    ct[0] = dom.nz_glob;
    ncwrap( ncmpi_begin_indep_data(ncid) , __LINE__ );
    ncwrap( ncmpi_put_vara_float( ncid , zVar    , st , ct , zCoord.data() ) , __LINE__ );
    // for (int k=0; k<dom.nz; k++) {
    Kokkos::parallel_for( dom.nz , KOKKOS_LAMBDA (int const k) {
      zCoord(k) = state.hyDensCells     (hs+k);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float( ncid , hyrVar  , st , ct , zCoord.data() ) , __LINE__ );
    // for (int k=0; k<dom.nz; k++) {
    Kokkos::parallel_for( dom.nz , KOKKOS_LAMBDA (int const k) {
      zCoord(k) = state.hyDensThetaCells(hs+k);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float( ncid , hyrtVar , st , ct , zCoord.data() ) , __LINE__ );
    ncwrap( ncmpi_end_indep_data(ncid) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void output(State &state, Domain const &dom, Parallel const &par) {
    int dimids[4];
    MPI_Offset st[1], ct[1];

    outTimer += dom.dt;
    if (outTimer < outFreq) {
      return;
    } else {
      outTimer -= outFreq;
    }

    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , "output.nc" , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "density" , &rVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "u"       , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "v"       , &vVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "w"       , &wVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "theta"   , &thVar ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "pressure", &pVar  ) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void writeState(State &state, Domain const &dom, Parallel const &par) {
    real3d data("data",dom.nz,dom.ny,dom.nx);
    MPI_Offset st[4], ct[4];

    st[0] = numOut; st[1] = 0     ; st[2] = par.j_beg; st[3] = par.i_beg;
    ct[0] = 1     ; ct[1] = dom.nz; ct[2] = dom.ny   ; ct[3] = dom.nx   ;

    // Write out density perturbation
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = state.state(idR,hs+k,hs+j,hs+i);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , rVar , st , ct , data.data() ) , __LINE__ );

    // Write out u wind
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = state.state(idRU,hs+k,hs+j,hs+i) / ( state.state(idR,hs+k,hs+j,hs+i) + state.hyDensCells(hs+k) );
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , uVar , st , ct , data.data() ) , __LINE__ );

    // Write out v wind
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = state.state(idRV,hs+k,hs+j,hs+i) / ( state.state(idR,hs+k,hs+j,hs+i) + state.hyDensCells(hs+k) );
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , vVar , st , ct , data.data() ) , __LINE__ );

    // Write out w wind
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = state.state(idRW,hs+k,hs+j,hs+i) / ( state.state(idR,hs+k,hs+j,hs+i) + state.hyDensCells(hs+k) );
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , wVar , st , ct , data.data() ) , __LINE__ );

    // Write out potential temperature perturbations
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = ( state.state(idRT,hs+k,hs+j,hs+i) + state.hyDensThetaCells(hs+k) ) /
                    ( state.state(idR ,hs+k,hs+j,hs+i) + state.hyDensCells     (hs+k) ) -
                    state.hyDensThetaCells(hs+k) / state.hyDensCells(hs+k);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , thVar , st , ct , data.data() ) , __LINE__ );

    // Write out perturbation pressure
    // for (int k=0; k<dom.nz; k++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nz*dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob) {
      int k, j, i;
      unpackIndices(iGlob,dom.nz,dom.ny,dom.nx,k,j,i);
      data(k,j,i) = C0*mypow(state.state(idRT,hs+k,hs+j,hs+i)+state.hyDensThetaCells(hs+k),GAMMA) -
                    C0*mypow(state.hyDensThetaCells(hs+k),GAMMA);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , pVar , st , ct , data.data() ) , __LINE__ );
  }


  //Error reporting routine for the PNetCDF I/O
  void ncwrap( int ierr , int line ) {
    if (ierr != NC_NOERR) {
      printf("NetCDF Error at line: %d\n", line);
      printf("%s\n",ncmpi_strerror(ierr));
      exit(-1);
    }
  }

};

#endif
