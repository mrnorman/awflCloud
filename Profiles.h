
#pragma once

#include "phys_params.h"

namespace profiles {

  YAKL_INLINE real initConstTheta_density(real const t0, real const z) {
    real exner = 1._fp - GRAV*z/(CP*t0);
    real p = pow( exner , CP/RD ) * P0;
    real rt = pow( p/C0 , 1._fp/GAMMA );
    return rt / t0;
  }


  YAKL_INLINE real initConstTheta_pressure(real const t0, real const z) {
    real r = initConstTheta_density(t0,z);
    return C0*pow(r*t0,GAMMA);
  }


  YAKL_INLINE real initConstTheta_pressureDeriv(real const t0, real const z) {
    real p = initConstTheta_pressure(t0,z);
    return -GRAV/(t0*RD)*pow(P0,RD/CP)*pow(p,-RD/CP+1);
  }


  YAKL_INLINE real initConstBVF_pot_temp(real const t0, real const bvf, real const z) {
    return t0 * exp(bvf*bvf*z/GRAV);
  }


  YAKL_INLINE real initConstBVF_density(real const t0, real const bvf, real const z) {
    real t = initConstBVF_pot_temp(t0,bvf,z);
    real exner = 1._fp - GRAV*GRAV/(CP*bvf*bvf) * (t-t0)/(t*t0);
    real p = pow( exner , CP/RD ) * P0;
    real rt = pow( p/C0 , 1._fp/GAMMA );
    return rt / t;
  }


  YAKL_INLINE real initConstBVF_pressure(real const t0, real const bvf, real const z) {
    real t = initConstBVF_pot_temp(t0,bvf,z);
    real r = initConstBVF_density (t0,bvf,z);
    return C0*pow(r*t,GAMMA);
  }


  YAKL_INLINE real ellipsoid_linear(real const x   , real const y   , real const z ,
                                    real const x0  , real const y0  , real const z0,
                                    real const xrad, real const yrad, real const zrad, real const amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    return amp * max( 1._fp - dist , 0._fp );
  }


  YAKL_INLINE real ellipsoid_cosine(real const x   , real const y   , real const z ,
                                    real const x0  , real const y0  , real const z0,
                                    real const xrad, real const yrad, real const zrad, real const amp, real const pwr) {
    real val = 0;
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    if (dist <= 1._fp) {
      val = amp * pow( (cos(M_PI*dist)+1)/2 , pwr );
    }
    return val;
  }


  YAKL_INLINE real igw(real const x   , real const y   , real const z ,
                       real const x0  , real const y0  , real const hr,
                       real const zlen, real const amp) {
    real xn = (x-x0)/hr;
    real yn = (y-y0)/hr;
    return amp * sin(M_PI*z/zlen) / sqrt( 1 + xn*xn + yn*yn );
  }

}

