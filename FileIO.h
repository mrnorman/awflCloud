
#ifndef _FILEIO_H_
#define _FILEIO_H_

#include "const.h"
#include "Array.h"
#include "State.h"
#include "pnetcdf.h"
#include "mpi.h"

class FileIO {

protected:

  int ncid, numOut;
  int tDim, xDim, yDim, zDim;
  int tVar, xVar, yVar, zVar, rVar, uVar, vVar, wVar, thVar, hyrVar, hyrtVar;

public:

  void outputInit(State &state, Domain const &dom, Parallel const &par) {
    int dimids[4];
    MPI_Offset st[1], ct[1];
    Array<real> xCoord(dom.nx);
    Array<real> yCoord(dom.ny);
    Array<real> zCoord(dom.nz);

    numOut = 0;

    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "z" , (MPI_Offset) dom.nz_glob  , &zDim ) , __LINE__ );

    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids , &tVar ) , __LINE__ );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_DOUBLE , 1 , dimids , &xVar ) , __LINE__ );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_DOUBLE , 1 , dimids , &yVar ) , __LINE__ );
    dimids[0] = zDim;
    ncwrap( ncmpi_def_var( ncid , "z"      , NC_DOUBLE , 1 , dimids , &zVar ) , __LINE__ );
    dimids[0] = tDim; dimids[1] = zDim; dimids[2] = yDim; dimids[3] = xDim;
    ncwrap( ncmpi_def_var( ncid , "density" , NC_DOUBLE , 4 , dimids , &rVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "u"       , NC_DOUBLE , 4 , dimids , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "v"       , NC_DOUBLE , 4 , dimids , &vVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "w"       , NC_DOUBLE , 4 , dimids , &wVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "theta"   , NC_DOUBLE , 4 , dimids , &thVar ) , __LINE__ );
    dimids[0] = zDim;
    ncwrap( ncmpi_def_var( ncid , "hyDens"      , NC_DOUBLE , 1 , dimids , &hyrVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "hyDensTheta" , NC_DOUBLE , 1 , dimids , &hyrtVar ) , __LINE__ );

    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Compute x, y, and z coordinates
    for (int i=0; i<dom.nx; i++) { xCoord(i) = ( par.i_beg + i + 0.5_fp ) * dom.dx; }
    for (int j=0; j<dom.ny; j++) { yCoord(j) = ( par.j_beg + j + 0.5_fp ) * dom.dy; }
    for (int k=0; k<dom.nz; k++) { zCoord(k) = (             k + 0.5_fp ) * dom.dz; }

    // Write out x, y, and z coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap( ncmpi_put_vara_double_all( ncid , xVar , st , ct , xCoord.get_data() ) , __LINE__ );
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap( ncmpi_put_vara_double_all( ncid , yVar , st , ct , yCoord.get_data() ) , __LINE__ );

    // Write out the hydrostatic background states and z coordinates
    st[0] = 0;
    ct[0] = dom.nz_glob;
    ncwrap( ncmpi_begin_indep_data(ncid) , __LINE__ );
    ncwrap( ncmpi_put_vara_double( ncid , zVar    , st , ct , zCoord                .get_data() ) , __LINE__ );
    ncwrap( ncmpi_put_vara_double( ncid , hyrVar  , st , ct , state.hyDensCells     .get_data() ) , __LINE__ );
    ncwrap( ncmpi_put_vara_double( ncid , hyrtVar , st , ct , state.hyDensThetaCells.get_data() ) , __LINE__ );
    ncwrap( ncmpi_end_indep_data(ncid) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void output(State &state, Domain const &dom, Parallel const &par) {
    int dimids[4];
    MPI_Offset st[1], ct[1];
    Array<real> xCoord(dom.nx);
    Array<real> yCoord(dom.ny);
    Array<real> zCoord(dom.nz);

    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , "output.nc" , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "density" , &rVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "u"       , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "v"       , &vVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "w"       , &wVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "theta"   , &thVar ) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void writeState(State &state, Domain const &dom, Parallel const &par) {
    Array<real> data(dom.nz,dom.ny,dom.nx);
    MPI_Offset st[4], ct[4];

    st[0] = numOut; st[1] = 0     ; st[2] = par.j_beg; st[3] = par.i_beg;
    ct[0] = 1     ; ct[1] = dom.nz; ct[2] = dom.ny   ; ct[3] = dom.nx   ;

    // Write out density perturbation
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          data(k,j,i) = state.state(idR,hs+k,hs+j,hs+i) - state.hyDensCells(k);
        }
      }
    }
    ncwrap( ncmpi_put_vara_double_all( ncid , rVar , st , ct , data.get_data() ) , __LINE__ );

    // Write out u wind
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          data(k,j,i) = state.state(idRU,hs+k,hs+j,hs+i) / state.state(idR,hs+k,hs+j,hs+i);
        }
      }
    }
    ncwrap( ncmpi_put_vara_double_all( ncid , uVar , st , ct , data.get_data() ) , __LINE__ );

    // Write out v wind
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          data(k,j,i) = state.state(idRV,hs+k,hs+j,hs+i) / state.state(idR,hs+k,hs+j,hs+i);
        }
      }
    }
    ncwrap( ncmpi_put_vara_double_all( ncid , vVar , st , ct , data.get_data() ) , __LINE__ );

    // Write out w wind
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          data(k,j,i) = state.state(idRW,hs+k,hs+j,hs+i) / state.state(idR,hs+k,hs+j,hs+i);
        }
      }
    }
    ncwrap( ncmpi_put_vara_double_all( ncid , wVar , st , ct , data.get_data() ) , __LINE__ );

    // Write out potential temperature perturbations
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          data(k,j,i) = state.state(idRT,hs+k,hs+j,hs+i) / state.state(idR,hs+k,hs+j,hs+i) -
                        state.hyDensThetaCells(k) / state.hyDensCells(k);
        }
      }
    }
    ncwrap( ncmpi_put_vara_double_all( ncid , thVar , st , ct , data.get_data() ) , __LINE__ );
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
