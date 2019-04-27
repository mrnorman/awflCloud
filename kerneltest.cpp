
#include "const.h"
#include "Array.h"
#include <iostream>

#ifdef _TARGET_GPU
#define _HOSTDEV __host__ __device__
#define _GLOBAL __global__
#else
#define _HOSTDEV 
#define _GLOBAL 
#endif

template <class F, class ... Ts> void launchCPU( int nIter , F const &func , Ts&&... args) {
  for (int i=0; i<nIter; i++) {
    func(i,args...);
  }
}

#ifdef _TARGET_GPU
template <class F, class ... Ts> _GLOBAL void launchCUDA( int nIter , F const &func , Ts&&... args) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < nIter) {
    func(i,args...);
  }
}
#endif

template <class F, class ... Ts> void launch( int nIter , F const &func , Ts&&... args) {
#ifdef _TARGET_GPU
  launchCUDA <<< nIter/128+1 , 128 >>> (nIter,func,args...);
#else
  launchCPU(nIter,func,args...);
#endif
}


#ifdef _TARGET_GPU
void synchronizeCUDA() { cudaDeviceSynchronize(); }
#endif

void synchronizeCPU() { }

void synchronize() {
#ifdef _TARGET_GPU
  synchronizeCUDA();
#else
  synchronizeCPU();
#endif
}


int main() {
  Array<real> a, b, c;
  int n = 1024*1024;
  a.setup(n);
  b.setup(n);
  c.setup(n);
  
  a = 2;
  b = 3;

  launch( n , [] _HOSTDEV (int i, Array<real> &a, Array<real> &b, Array<real> &c) { c(i) = a(i) + b(i); } , a , b , c );
  synchronize();

  std::cout << (int) c.sum() << " " << 5*1024*1024 << "\n";
}

