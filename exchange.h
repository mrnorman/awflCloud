
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "mpi.h"

namespace exchange {

  int const maxPack = numState*2;

  MPI_Request sReq [2];
  MPI_Request rReq [2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  real4d haloSendBufS;
  real4d haloSendBufN;
  real4d haloSendBufW;
  real4d haloSendBufE;
  real4d haloRecvBufS;
  real4d haloRecvBufN;
  real4d haloRecvBufW;
  real4d haloRecvBufE;

  real3d edgeRecvBufE;
  real3d edgeRecvBufW;
  real3d edgeSendBufE;
  real3d edgeSendBufW;
  real3d edgeRecvBufN;
  real3d edgeRecvBufS;
  real3d edgeSendBufN;
  real3d edgeSendBufS;


  inline void allocate(Domain &dom) {
    haloSendBufS = real4d("haloSendBufS",maxPack,dom.nz,hs,dom.nx);
    haloSendBufN = real4d("haloSendBufN",maxPack,dom.nz,hs,dom.nx);
    haloSendBufW = real4d("haloSendBufW",maxPack,dom.nz,dom.ny,hs);
    haloSendBufE = real4d("haloSendBufE",maxPack,dom.nz,dom.ny,hs);
    haloRecvBufS = real4d("haloRecvBufS",maxPack,dom.nz,hs,dom.nx);
    haloRecvBufN = real4d("haloRecvBufN",maxPack,dom.nz,hs,dom.nx);
    haloRecvBufW = real4d("haloRecvBufW",maxPack,dom.nz,dom.ny,hs);
    haloRecvBufE = real4d("haloRecvBufE",maxPack,dom.nz,dom.ny,hs);

    edgeSendBufS = real3d("edgeSendBufS",maxPack,dom.nz,dom.nx);
    edgeSendBufN = real3d("edgeSendBufN",maxPack,dom.nz,dom.nx);
    edgeSendBufW = real3d("edgeSendBufW",maxPack,dom.nz,dom.ny);
    edgeSendBufE = real3d("edgeSendBufE",maxPack,dom.nz,dom.ny);
    edgeRecvBufS = real3d("edgeRecvBufS",maxPack,dom.nz,dom.nx);
    edgeRecvBufN = real3d("edgeRecvBufN",maxPack,dom.nz,dom.nx);
    edgeRecvBufW = real3d("edgeRecvBufW",maxPack,dom.nz,dom.ny);
    edgeRecvBufE = real3d("edgeRecvBufE",maxPack,dom.nz,dom.ny);
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }

  inline void haloPackN_x(Domain const &dom, real4d const &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    //       for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j, ii;
      unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
      haloSendBufW(nPack+v,k,j,ii) = a(v,hs+k,hs+j,hs    +ii);
      haloSendBufE(nPack+v,k,j,ii) = a(v,hs+k,hs+j,dom.nx+ii);
    });
    nPack = nPack + n;
  }


  inline void haloPackN_y(Domain const &dom, real4d const &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int ii=0; ii<hs; ii++) {
    //       for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, ii, i;
      unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
      haloSendBufS(nPack+v,k,ii,i) = a(v,hs+k,hs    +ii,hs+i);
      haloSendBufN(nPack+v,k,ii,i) = a(v,hs+k,dom.ny+ii,hs+i);
    });
    nPack = nPack + n;
  }


  // inline void haloPack1_x(Domain const &dom, real3d const &a) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int ii=0; ii<hs; ii++) {
  //         haloSendBufW(nPack,k,j,ii) = a(hs+k,hs+j,hs    +ii);
  //         haloSendBufE(nPack,k,j,ii) = a(hs+k,hs+j,dom.nx+ii);
  //       }
  //     }
  //   }
  //   nPack = nPack + 1;
  // }


  // inline void haloPack1_y(Domain const &dom, real3d const &a) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int ii=0; ii<hs; ii++) {
  //       for (int i=0; i<dom.nx; i++) {
  //         haloSendBufS(nPack,k,ii,i) = a(hs+k,hs    +ii,hs+i);
  //         haloSendBufN(nPack,k,ii,i) = a(hs+k,dom.ny+ii,hs+i);
  //       }
  //     }
  //   }
  //   nPack = nPack + 1;
  // }


  inline void haloUnpackN_x(Domain const &dom, real4d &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    //       for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j, ii;
      unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
      a(v,hs+k,hs+j,          ii) = haloRecvBufW(nUnpack+v,k,j,ii);
      a(v,hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack+v,k,j,ii);
    });
    nUnpack = nUnpack + n;
  }


  inline void haloUnpackN_y(Domain const &dom, real4d &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int ii=0; ii<hs; ii++) {
    //       for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, ii, i;
      unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
      a(v,hs+k,          ii,hs+i) = haloRecvBufS(nUnpack+v,k,ii,i);
      a(v,hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack+v,k,ii,i);
    });
    nUnpack = nUnpack + n;
  }


  // inline void haloUnpack1_x(Domain const &dom, real3d &a) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  //       for (int ii=0; ii<hs; ii++) {
  //         a(hs+k,hs+j,          ii) = haloRecvBufW(nUnpack,k,j,ii);
  //         a(hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack,k,j,ii);
  //       }
  //     }
  //   }
  //   nUnpack = nUnpack + 1;
  // }


  // inline void haloUnpack1_y(Domain const &dom, real3d &a) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int ii=0; ii<hs; ii++) {
  //       for (int i=0; i<dom.nx; i++) {
  //         a(hs+k,          ii,hs+i) = haloRecvBufS(nUnpack,k,ii,i);
  //         a(hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack,k,ii,i);
  //       }
  //     }
  //   }
  //   nUnpack = nUnpack + 1;
  // }


  inline void haloExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    Kokkos::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufW.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufE.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufW.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufE.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void haloExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    Kokkos::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufS.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufN.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufS.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufN.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgePackN_x(Domain const &dom, real5d const &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j;
      unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
      edgeSendBufW(nPack+v,k,j) = a(v,1,k,j,0     );
      edgeSendBufE(nPack+v,k,j) = a(v,0,k,j,dom.nx);
    });
    nPack = nPack + n;
  }


  inline void edgePackN_y(Domain const &dom, real5d const &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, i;
      unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
      edgeSendBufS(nPack+v,k,i) = a(v,1,k,0     ,i);
      edgeSendBufN(nPack+v,k,i) = a(v,0,k,dom.ny,i);
    });
    nPack = nPack + n;
  }


  inline void edgeUnpackN_x(Domain const &dom, real5d &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j;
      unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
      a(v,0,k,j,0     ) = edgeRecvBufW(nUnpack+v,k,j);
      a(v,1,k,j,dom.nx) = edgeRecvBufE(nUnpack+v,k,j);
    });
    nUnpack = nUnpack + n;
  }


  inline void edgeUnpackN_y(Domain const &dom, real5d &a, int const n) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, i;
      unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
      a(v,0,k,0     ,i) = edgeRecvBufS(nUnpack+v,k,i);
      a(v,1,k,dom.ny,i) = edgeRecvBufN(nUnpack+v,k,i);
    });
    nUnpack = nUnpack + n;
  }


  inline void edgeExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    Kokkos::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufW.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufE.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufW.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufE.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgeExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    Kokkos::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufS.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufN.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufS.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufN.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }

}

#endif
