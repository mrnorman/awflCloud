
#ifndef _CONST_H_
#define _CONST_H_

#include <cmath>
#include "YAKL.h"

using yakl::SArray;
using yakl::Array;
using yakl::styleC;
using yakl::memHost;
using yakl::memDevice;
using yakl::c::Bounds;
using yakl::c::parallel_for;

typedef float         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

typedef Array<real,1,memDevice,styleC> real1d;
typedef Array<real,2,memDevice,styleC> real2d;
typedef Array<real,3,memDevice,styleC> real3d;
typedef Array<real,4,memDevice,styleC> real4d;
typedef Array<real,5,memDevice,styleC> real5d;
typedef Array<real,6,memDevice,styleC> real6d;
typedef Array<real,7,memDevice,styleC> real7d;
typedef Array<real,8,memDevice,styleC> real8d;

typedef Array<real,1,memHost,styleC> realHost1d;
typedef Array<real,2,memHost,styleC> realHost2d;
typedef Array<real,3,memHost,styleC> realHost3d;
typedef Array<real,4,memHost,styleC> realHost4d;
typedef Array<real,5,memHost,styleC> realHost5d;
typedef Array<real,6,memHost,styleC> realHost6d;
typedef Array<real,7,memHost,styleC> realHost7d;
typedef Array<real,8,memHost,styleC> realHost8d;

#include "params.h"

YAKL_INLINE real constexpr operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

int constexpr ord      = 9;
int constexpr tord     = 5;
int constexpr hs       = (ord-1)/2;
int constexpr numState = 5;

int constexpr idR      = 0;
int constexpr idRU     = 1;
int constexpr idRV     = 2;
int constexpr idRW     = 3;
int constexpr idRT     = 4;

int constexpr idU      = 1;
int constexpr idV      = 2;
int constexpr idW      = 3;
int constexpr idT      = 4;

// Some physical constants
real constexpr PI    = 3.1415926535897932384626433832795028842;
real constexpr GRAV  = 9.8;
real constexpr CP    = 1004.;
real constexpr CV    = 717.;
real constexpr RD    = 287.;
real constexpr P0    = 1.0e5;
real constexpr C0    = 27.5629410929725921310572974482;
real constexpr GAMMA  = 1.40027894002789400278940027894;

template <class T> YAKL_INLINE T mymin( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> YAKL_INLINE T mymax( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

#endif
