
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "const.h"

class Domain {

public:

  ulong nx_glob;
  ulong ny_glob;
  ulong nz_glob;

  uint nx;
  uint ny;
  uint nz;

  int run2d;

  real xlen;
  real ylen;
  real zlen;

  real dx;
  real dy;
  real dz;

  real cfl;
  real simLength;

  real etime;

  real dt;
};

#endif
