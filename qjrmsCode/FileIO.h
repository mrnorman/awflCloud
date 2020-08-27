
#ifndef _FILEIO_H_
#define _FILEIO_H_

#include "const.h"
#include "pnetcdf.h"
#include "mpi.h"
#include "YAKL.h"
#include "Domain.h"
#include "Parallel.h"

class FileIO {

protected:

  real outTimer;
  int ncid, numOut;
  int tDim, xDim, yDim, zDim;
  int tVar, xVar, yVar, zVar, rVar, uVar, vVar, wVar, thVar, hyrVar, hyrtVar, hypVar, pVar;

public:

  void outputInit(real4d &state, Domain const &dom, Parallel const &par);

  void output(real4d &state, Domain const &dom, Parallel const &par);

  void writeState(real4d &state, Domain const &dom, Parallel const &par);

  void writeStateThetaPrim(real4d &state, Domain const &dom, Parallel const &par);

  void writeStateThetaCons(real4d &state, Domain const &dom, Parallel const &par);

  //Error reporting routine for the PNetCDF I/O
  void ncwrap( int ierr , int line );
};

#endif
