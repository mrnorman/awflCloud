
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "Array.h"
#include "mpi.h"

class Exchange {

protected:

  int const maxPack = numState*2;

  MPI_Request sReq [2];
  MPI_Request rReq [2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  Array<real> haloSendBufS;
  Array<real> haloSendBufN;
  Array<real> haloSendBufW;
  Array<real> haloSendBufE;
  Array<real> haloRecvBufS;
  Array<real> haloRecvBufN;
  Array<real> haloRecvBufW;
  Array<real> haloRecvBufE;

  Array<real> edgeRecvBufE;
  Array<real> edgeRecvBufW;
  Array<real> edgeSendBufE;
  Array<real> edgeSendBufW;
  Array<real> edgeRecvBufN;
  Array<real> edgeRecvBufS;
  Array<real> edgeSendBufN;
  Array<real> edgeSendBufS;

public:


  inline void initialize(Domain &dom) {
    haloSendBufS.setup(maxPack,dom.nz,hs,dom.nx);
    haloSendBufN.setup(maxPack,dom.nz,hs,dom.nx);
    haloSendBufW.setup(maxPack,dom.nz,dom.ny,hs);
    haloSendBufE.setup(maxPack,dom.nz,dom.ny,hs);
    haloRecvBufS.setup(maxPack,dom.nz,hs,dom.nx);
    haloRecvBufN.setup(maxPack,dom.nz,hs,dom.nx);
    haloRecvBufW.setup(maxPack,dom.nz,dom.ny,hs);
    haloRecvBufE.setup(maxPack,dom.nz,dom.ny,hs);

    edgeSendBufS.setup(maxPack,dom.nz,dom.nx);
    edgeSendBufN.setup(maxPack,dom.nz,dom.nx);
    edgeSendBufW.setup(maxPack,dom.nz,dom.ny);
    edgeSendBufE.setup(maxPack,dom.nz,dom.ny);
    edgeRecvBufS.setup(maxPack,dom.nz,dom.nx);
    edgeRecvBufN.setup(maxPack,dom.nz,dom.nx);
    edgeRecvBufW.setup(maxPack,dom.nz,dom.ny);
    edgeRecvBufE.setup(maxPack,dom.nz,dom.ny);
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }

  inline void haloPackN_x(Domain const &dom, Array<real> const &a, int const n) {
    launcher.parallelFor( n*dom.nz*dom.ny*hs , 
      [this] _YAKL (int iGlob, Domain const &dom, Array<real> const &a, int const n) {
        int v, k, j, ii;
        yakl::unpackIndices(iGlob, n, dom.nz, dom.ny, hs, v, k, j, ii);
        haloSendBufW(nPack+v,k,j,ii) = a(v,hs+k,hs+j,hs    +ii);
        haloSendBufE(nPack+v,k,j,ii) = a(v,hs+k,hs+j,dom.nx+ii);
      } , dom , a , n );
    launcher.synchronizeSelf();
  }


  inline void haloPackN_y(Domain const &dom, Array<real> const &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int ii=0; ii<hs; ii++) {
          for (int i=0; i<dom.nx; i++) {
            haloSendBufS(nPack+v,k,ii,i) = a(v,hs+k,hs    +ii,hs+i);
            haloSendBufN(nPack+v,k,ii,i) = a(v,hs+k,dom.ny+ii,hs+i);
          }
        }
      }
    }
    nPack = nPack + n;
  }


  inline void haloPack1_x(Domain const &dom, Array<real> const &a) {
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int ii=0; ii<hs; ii++) {
          haloSendBufW(nPack,k,j,ii) = a(hs+k,hs+j,hs    +ii);
          haloSendBufE(nPack,k,j,ii) = a(hs+k,hs+j,dom.nx+ii);
        }
      }
    }
    nPack = nPack + 1;
  }


  inline void haloPack1_y(Domain const &dom, Array<real> const &a) {
    for (int k=0; k<dom.nz; k++) {
      for (int ii=0; ii<hs; ii++) {
        for (int i=0; i<dom.nx; i++) {
          haloSendBufS(nPack,k,ii,i) = a(hs+k,hs    +ii,hs+i);
          haloSendBufN(nPack,k,ii,i) = a(hs+k,dom.ny+ii,hs+i);
        }
      }
    }
    nPack = nPack + 1;
  }


  inline void haloUnpackN_x(Domain const &dom, Array<real> &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          for (int ii=0; ii<hs; ii++) {
            a(v,hs+k,hs+j,          ii) = haloRecvBufW(nUnpack+v,k,j,ii);
            a(v,hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack+v,k,j,ii);
          }
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void haloUnpackN_y(Domain const &dom, Array<real> &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int ii=0; ii<hs; ii++) {
          for (int i=0; i<dom.nx; i++) {
            a(v,hs+k,          ii,hs+i) = haloRecvBufS(nUnpack+v,k,ii,i);
            a(v,hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack+v,k,ii,i);
          }
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void haloUnpack1_x(Domain const &dom, Array<real> &a) {
    for (int k=0; k<dom.nz; k++) {
      for (int j=0; j<dom.ny; j++) {
        for (int ii=0; ii<hs; ii++) {
          a(hs+k,hs+j,          ii) = haloRecvBufW(nUnpack,k,j,ii);
          a(hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack,k,j,ii);
        }
      }
    }
    nUnpack = nUnpack + 1;
  }


  inline void haloUnpack1_y(Domain const &dom, Array<real> &a) {
    for (int k=0; k<dom.nz; k++) {
      for (int ii=0; ii<hs; ii++) {
        for (int i=0; i<dom.nx; i++) {
          a(hs+k,          ii,hs+i) = haloRecvBufS(nUnpack,k,ii,i);
          a(hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack,k,ii,i);
        }
      }
    }
    nUnpack = nUnpack + 1;
  }


  inline void haloExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufW.get_data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufE.get_data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufW.get_data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufE.get_data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void haloExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufS.get_data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufN.get_data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufS.get_data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufN.get_data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgePackN_x(Domain const &dom, Array<real> const &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          edgeSendBufW(nPack+v,k,j) = a(v,1,k,j,0     );
          edgeSendBufE(nPack+v,k,j) = a(v,0,k,j,dom.nx);
        }
      }
    }
    nPack = nPack + n;
  }


  inline void edgePackN_y(Domain const &dom, Array<real> const &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int i=0; i<dom.nx; i++) {
          edgeSendBufS(nPack+v,k,i) = a(v,1,k,0     ,i);
          edgeSendBufN(nPack+v,k,i) = a(v,0,k,dom.ny,i);
        }
      }
    }
    nPack = nPack + n;
  }


  inline void edgeUnpackN_x(Domain const &dom, Array<real> &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int j=0; j<dom.ny; j++) {
          a(v,0,k,j,0     ) = edgeRecvBufW(nUnpack+v,k,j);
          a(v,1,k,j,dom.nx) = edgeRecvBufE(nUnpack+v,k,j);
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void edgeUnpackN_y(Domain const &dom, Array<real> &a, int const n) {
    for (int v=0; v<n; v++) {
      for (int k=0; k<dom.nz; k++) {
        for (int i=0; i<dom.nx; i++) {
          a(v,0,k,0     ,i) = edgeRecvBufS(nUnpack+v,k,i);
          a(v,1,k,dom.ny,i) = edgeRecvBufN(nUnpack+v,k,i);
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void edgeExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufW.get_data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufE.get_data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufW.get_data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufE.get_data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgeExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufS.get_data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufN.get_data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufS.get_data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufN.get_data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }

};

#endif
