
#include "const.h"
#include "Array.h"
#include <iostream>
#include <type_traits>

class Nope {
};
Nope _nope;

template <class F, class T1=Nope, class T2=Nope, class T3=Nope, class T4=Nope, class T5=Nope, 
                   class T6=Nope, class T7=Nope, class T8=Nope, class T9=Nope, class T10=Nope>
  void launch( ulong nIter , F const &func , T1 &p1=_nope , T2 &p2=_nope , T3 &p3=_nope , T4 &p4=_nope , T5  &p5 =_nope , 
                                             T6 &p6=_nope , T7 &p7=_nope , T8 &p8=_nope , T9 &p9=_nope , T10 &p10=_nope ) {
  for (ulong i=0; i<nIter; i++) {
    if        constexpr ( std::is_same<T1 ,Nope>::value ) {
      func(i);
    } else if constexpr ( std::is_same<T2 ,Nope>::value ) {
      func(i,p1);
    } else if constexpr ( std::is_same<T3 ,Nope>::value ) {
      func(i,p1,p2);
    } else if constexpr ( std::is_same<T4 ,Nope>::value ) {
      func(i,p1,p2,p3);
    } else if constexpr ( std::is_same<T5 ,Nope>::value ) {
      func(i,p1,p2,p3,p4);
    } else if constexpr ( std::is_same<T6 ,Nope>::value ) {
      func(i,p1,p2,p3,p4,p5);
    } else if constexpr ( std::is_same<T7 ,Nope>::value ) {
      func(i,p1,p2,p3,p4,p5,p6);
    } else if constexpr ( std::is_same<T8 ,Nope>::value ) {
      func(i,p1,p2,p3,p4,p5,p6,p7);
    } else if constexpr ( std::is_same<T9 ,Nope>::value ) {
      func(i,p1,p2,p3,p4,p5,p6,p7,p8);
    } else if constexpr ( std::is_same<T10,Nope>::value ) {
      func(i,p1,p2,p3,p4,p5,p6,p7,p8,p9);
    } else {
      func(i,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
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
