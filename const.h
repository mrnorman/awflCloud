
#ifndef _CONST_H_
#define _CONST_H_

#include "math.h"

typedef float         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

#include "params.h"

inline _YAKL real operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

#define ord  5
#define tord 3
int const hs = (ord-1)/2;

int const numState = 5;

int const idR  = 0;
int const idRU = 1;
int const idRV = 2;
int const idRW = 3;
int const idRT = 4;

// Some physical constants
real const PI    = 3.1415926535897932384626433832795028842;
real const GRAV  = 9.8;
real const CP    = 1004.;
real const CV    = 717.;
real const RD    = 287.;
real const P0    = 1.0e5;
real const C0    = 27.5629410929725921310572974482;
real const GAMMA  = 1.40027894002789400278940027894;

inline _YAKL double mypow ( double const x , double const p ) { return pow (x,p); }
inline _YAKL float  mypow ( float  const x , float  const p ) { return powf(x,p); }
inline _YAKL double mysqrt( double const x ) { return sqrt (x); }
inline _YAKL float  mysqrt( float  const x ) { return sqrtf(x); }
inline _YAKL double myfabs( double const x ) { return fabs (x); }
inline _YAKL float  myfabs( float  const x ) { return fabsf(x); }

template <class T> inline _YAKL T min( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> inline _YAKL T max( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

#endif
