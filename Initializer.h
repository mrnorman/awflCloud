
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Hydrostasis.h"
#include "Exchange.h"
#include "TransformMatrices.h"
#include "TimeIntegrator.h"
#include "mpi.h"
#include "Domain.h"
#include "Parallel.h"

class Random {
public:
  uint64_t state;

  // No seed (1.5 billion-th prime number)
  YAKL_INLINE Random()                            { this->state = (uint64_t) 3121238909U;}
  YAKL_INLINE Random(uint64_t state)              { this->state = state; warmup(); }
  YAKL_INLINE Random(Random const            &in) { this->state = in.state; }
  YAKL_INLINE Random(Random                 &&in) { this->state = in.state; }
  YAKL_INLINE Random &operator=(Random const &in) { this->state = in.state; return *this; }
  YAKL_INLINE Random &operator=(Random      &&in) { this->state = in.state; return *this; }

  // Warmup (probably not needed, but whatever)
  YAKL_INLINE void warmup(int ncycles=10) {
    for (int i=0; i<ncycles; i++) { gen(); }
  }

  // Return a random unsigned 32-bit integer
  YAKL_INLINE uint32_t gen() {
      uint32_t c = state>>32;
      uint32_t x = state&0xFFFFFFFF;
      state = x*((uint64_t)4294883355U) + c;
      return x^c;
  }

  // Return floating point value: domain \in (0,1]
  template <class T> YAKL_INLINE T genFP() {
    return ( (T) gen() + 1 ) / ( (T) 4294967295U );
  }

  // Return a floating point value with custom bounds
  template <class T> YAKL_INLINE T genFP(T lb, T ub) {
    return  genFP<T>() * (ub-lb) + lb;
  }

  // Return floating point value: domain \in (0,1]
  template <class T> YAKL_INLINE void fillArray(T *data, int n) {
    for (int i=0; i<n; i++) {
      data[i] = genFP<T>();
    }
  }

};




void initialize_mpi( int *argc , char ***argv , Parallel &par );

void initialize(real4d &state, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint);

YAKL_INLINE real ellipsoid_linear(real const x   , real const y   , real const z ,
                                  real const x0  , real const y0  , real const z0,
                                  real const xrad, real const yrad, real const zrad, real const amp) {
  real xn = (x-x0)/xrad;
  real yn = (y-y0)/yrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + yn*yn + zn*zn );
  return amp * mymax( 1._fp - dist , 0._fp );
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
    val = amp * pow( (cos(PI*dist)+1)/2 , pwr );
  }
  return val;
}


#endif
