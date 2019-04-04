
#ifndef _HYDROSTASIS_H_
#define _HYDROSTASIS_H_

#include "const.h"

class Hydrostasis {

public:

  /*
    Computes the hydrostatic density (r) at the given constant
    potential temperature (t0) and location (zloc)
  */
  inline void hydroConstTheta(real const t0, real const z, real &r) {
    real exner = 1._fp - GRAV*z/(CP*t0);
    real p = mypow( exner , CP/RD ) * P0;
    real rt = mypow( p/C0 , 1._fp/GAMMA );
    r = rt / t0;
  }

};

#endif
