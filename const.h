
#ifndef _CONST_H_
#define _CONST_H_

#ifdef __NVCC__
#define _HOSTDEV __host__ __device__
#else
#define _HOSTDEV
#endif

typedef float         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

inline _HOSTDEV real operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

#define ord  5
#define tord 2
uint const hs = (ord-1)/2;

#endif
