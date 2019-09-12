
#ifndef _CONST_H_
#define _CONST_H_

#include <cmath>
#include "Array.h"

typedef float         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

#if defined(__NVCC__) || defined(__USE_HIP__)
typedef yakl::Array<real,yakl::memDevice> realArr;
#else
typedef yakl::Array<real,yakl::memHost> realArr;
#endif

typedef yakl::Array<real,yakl::memHost> realArrHost;

#ifdef __NVCC__
#define _HOSTDEV __host__ __device__
#elif defined(__USE_HIP__)
#define _HOSTDEV __host__ __device__
#include "hip/hip_runtime.h"
#else
#define _HOSTDEV 
#endif

#include "params.h"

inline _HOSTDEV real operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

#define ord  9
#define tord 3
#define hs (ord-1)/2

#define numState 5

#define idR  0
#define idRU 1
#define idRV 2
#define idRW 3
#define idRT 4

#define idU 1
#define idV 2
#define idW 3
#define idT 4

// Some physical constants
real const PI    = 3.1415926535897932384626433832795028842;
real const GRAV  = 9.8;
real const CP    = 1004.;
real const CV    = 717.;
real const RD    = 287.;
real const P0    = 1.0e5;
real const C0    = 27.5629410929725921310572974482;
real const GAMMA  = 1.40027894002789400278940027894;

template <class T> inline _HOSTDEV T mymin( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> inline _HOSTDEV T mymax( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

#endif
