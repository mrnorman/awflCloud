
#include "const.h"
#include "Array.h"
#include <iostream>
#include <type_traits>

class noParam {
};
noParam _nope;

template <class F, class T1=noParam, class T2=noParam, class T3=noParam, class T4=noParam, class T5=noParam, class T6=noParam>
  void launch( ulong nIter , F const &func , T1 &p1=_nope , T2 &p2=_nope , T3 &p3=_nope , T4 &p4=_nope, T5 &p5=_nope, T6 &p6=_nope ) {
  for (ulong i=0; i<nIter; i++) {
    if        constexpr ( std::is_same<T1,noParam>::value ) {
      func(i);
    } else if constexpr ( std::is_same<T2,noParam>::value ) {
      func(i,p1);
    } else if constexpr ( std::is_same<T3,noParam>::value ) {
      func(i,p1,p2);
    } else if constexpr ( std::is_same<T4,noParam>::value ) {
      func(i,p1,p2,p3);
    } else if constexpr ( std::is_same<T5,noParam>::value ) {
      func(i,p1,p2,p3,p4);
    } else if constexpr ( std::is_same<T6,noParam>::value ) {
      func(i,p1,p2,p3,p4,p5);
    } else {
      func(i,p1,p2,p3,p4,p5,p6);
    }
  }
}

int main() {
  Array<float> a, b, c;
  ulong n = 1024*1024;
  a.setup(n);
  b.setup(n);
  c.setup(n);
  
  a = 2;
  b = 3;

  launch( n , []( ulong i , Array<float> &a , Array<float> &b , Array<float> &c ) { c(i) = a(i) + b(i); } , a , b , c );

  std::cout << (int) c.sum() << " " << 5*1024*1024 << "\n";
}
