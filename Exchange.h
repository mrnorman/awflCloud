
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
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

  real1d haloSendBufS;
  real1d haloSendBufN;
  real1d haloSendBufW;
  real1d haloSendBufE;
  real1d haloRecvBufS;
  real1d haloRecvBufN;
  real1d haloRecvBufW;
  real1d haloRecvBufE;

  real1d edgeRecvBufE;
  real1d edgeRecvBufW;
  real1d edgeSendBufE;
  real1d edgeSendBufW;
  real1d edgeRecvBufN;
  real1d edgeRecvBufS;
  real1d edgeSendBufN;
  real1d edgeSendBufS;

  real *haloSendBufS_cpu;
  real *haloSendBufN_cpu;
  real *haloSendBufW_cpu;
  real *haloSendBufE_cpu;
  real *haloRecvBufS_cpu;
  real *haloRecvBufN_cpu;
  real *haloRecvBufW_cpu;
  real *haloRecvBufE_cpu;

  real *edgeRecvBufE_cpu;
  real *edgeRecvBufW_cpu;
  real *edgeSendBufE_cpu;
  real *edgeSendBufW_cpu;
  real *edgeRecvBufN_cpu;
  real *edgeRecvBufS_cpu;
  real *edgeSendBufN_cpu;
  real *edgeSendBufS_cpu;

public:


  inline void allocate(Domain &dom) {
    haloSendBufS = real1d("haloSendBufS",maxPack*dom.nz*hs*dom.nx);
    haloSendBufN = real1d("haloSendBufN",maxPack*dom.nz*hs*dom.nx);
    haloSendBufW = real1d("haloSendBufW",maxPack*dom.nz*dom.ny*hs);
    haloSendBufE = real1d("haloSendBufE",maxPack*dom.nz*dom.ny*hs);
    haloRecvBufS = real1d("haloRecvBufS",maxPack*dom.nz*hs*dom.nx);
    haloRecvBufN = real1d("haloRecvBufN",maxPack*dom.nz*hs*dom.nx);
    haloRecvBufW = real1d("haloRecvBufW",maxPack*dom.nz*dom.ny*hs);
    haloRecvBufE = real1d("haloRecvBufE",maxPack*dom.nz*dom.ny*hs);

    edgeSendBufS = real1d("edgeSendBufS",maxPack*dom.nz*dom.nx);
    edgeSendBufN = real1d("edgeSendBufN",maxPack*dom.nz*dom.nx);
    edgeSendBufW = real1d("edgeSendBufW",maxPack*dom.nz*dom.ny);
    edgeSendBufE = real1d("edgeSendBufE",maxPack*dom.nz*dom.ny);
    edgeRecvBufS = real1d("edgeRecvBufS",maxPack*dom.nz*dom.nx);
    edgeRecvBufN = real1d("edgeRecvBufN",maxPack*dom.nz*dom.nx);
    edgeRecvBufW = real1d("edgeRecvBufW",maxPack*dom.nz*dom.ny);
    edgeRecvBufE = real1d("edgeRecvBufE",maxPack*dom.nz*dom.ny);

    #ifdef __NVCC__
      cudaMallocHost( &haloSendBufS_cpu , maxPack*dom.nz*hs*dom.nx*sizeof(real) );
      cudaMallocHost( &haloSendBufN_cpu , maxPack*dom.nz*hs*dom.nx*sizeof(real) );
      cudaMallocHost( &haloSendBufW_cpu , maxPack*dom.nz*dom.ny*hs*sizeof(real) );
      cudaMallocHost( &haloSendBufE_cpu , maxPack*dom.nz*dom.ny*hs*sizeof(real) );
      cudaMallocHost( &haloRecvBufS_cpu , maxPack*dom.nz*hs*dom.nx*sizeof(real) );
      cudaMallocHost( &haloRecvBufN_cpu , maxPack*dom.nz*hs*dom.nx*sizeof(real) );
      cudaMallocHost( &haloRecvBufW_cpu , maxPack*dom.nz*dom.ny*hs*sizeof(real) );
      cudaMallocHost( &haloRecvBufE_cpu , maxPack*dom.nz*dom.ny*hs*sizeof(real) );

      cudaMallocHost( &edgeSendBufS_cpu , maxPack*dom.nz*dom.nx*sizeof(real) );
      cudaMallocHost( &edgeSendBufN_cpu , maxPack*dom.nz*dom.nx*sizeof(real) );
      cudaMallocHost( &edgeSendBufW_cpu , maxPack*dom.nz*dom.ny*sizeof(real) );
      cudaMallocHost( &edgeSendBufE_cpu , maxPack*dom.nz*dom.ny*sizeof(real) );
      cudaMallocHost( &edgeRecvBufS_cpu , maxPack*dom.nz*dom.nx*sizeof(real) );
      cudaMallocHost( &edgeRecvBufN_cpu , maxPack*dom.nz*dom.nx*sizeof(real) );
      cudaMallocHost( &edgeRecvBufW_cpu , maxPack*dom.nz*dom.ny*sizeof(real) );
      cudaMallocHost( &edgeRecvBufE_cpu , maxPack*dom.nz*dom.ny*sizeof(real) );
    #else
      haloSendBufS_cpu = haloSendBufS.data(); 
      haloSendBufN_cpu = haloSendBufN.data(); 
      haloSendBufW_cpu = haloSendBufW.data(); 
      haloSendBufE_cpu = haloSendBufE.data(); 
      haloRecvBufS_cpu = haloRecvBufS.data(); 
      haloRecvBufN_cpu = haloRecvBufN.data(); 
      haloRecvBufW_cpu = haloRecvBufW.data(); 
      haloRecvBufE_cpu = haloRecvBufE.data(); 

      edgeSendBufS_cpu = edgeSendBufS.data(); 
      edgeSendBufN_cpu = edgeSendBufN.data(); 
      edgeSendBufW_cpu = edgeSendBufW.data(); 
      edgeSendBufE_cpu = edgeSendBufE.data(); 
      edgeRecvBufS_cpu = edgeRecvBufS.data(); 
      edgeRecvBufN_cpu = edgeRecvBufN.data(); 
      edgeRecvBufW_cpu = edgeRecvBufW.data(); 
      edgeRecvBufE_cpu = edgeRecvBufE.data(); 
    #endif
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }


  inline void haloPackN_x(Domain const &dom, real4d const &a, int const n) {
    haloPackN_x_ext(dom,a,n,haloSendBufW,haloSendBufE,nPack);
  }
  inline void haloPackN_x_ext(Domain const &dom, real4d const &a, int const n, real1d &haloSendBufW, real1d &haloSendBufE, int &nPack) {
    Kokkos::parallel_for( n*dom.nz*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j, ii;
      unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
      int nGlob = dom.nz*dom.ny*hs;
      haloSendBufW(nPack*nGlob+iGlob) = a(v,hs+k,hs+j,hs    +ii);
      haloSendBufE(nPack*nGlob+iGlob) = a(v,hs+k,hs+j,dom.nx+ii);
    });
    nPack = nPack + n;
  }


  inline void haloPackN_y(Domain const &dom, real4d const &a, int const n) {
    haloPackN_y_ext(dom,a,n,haloSendBufS,haloSendBufN,nPack);
  }
  inline void haloPackN_y_ext(Domain const &dom, real4d const &a, int const n, real1d &haloSendBufS, real1d &haloSendBufN, int &nPack) {
    Kokkos::parallel_for( n*dom.nz*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, ii, i;
      unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
      int nGlob = dom.nz*hs*dom.nx;
      haloSendBufS(nPack*nGlob+iGlob) = a(v,hs+k,hs    +ii,hs+i);
      haloSendBufN(nPack*nGlob+iGlob) = a(v,hs+k,dom.ny+ii,hs+i);
    });
    nPack = nPack + n;
  }


  inline void haloUnpackN_x(Domain const &dom, real4d &a, int const n) {
    haloUnpackN_x_ext(dom, a, n, haloRecvBufW, haloRecvBufE, nUnpack);
  }
  inline void haloUnpackN_x_ext(Domain const &dom, real4d &a, int const n, real1d const &haloRecvBufW, real1d const &haloRecvBufE, int &nUnpack) {
    Kokkos::parallel_for( n*dom.nz*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j, ii;
      unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
      int nGlob = dom.nz*dom.ny*hs;
      a(v,hs+k,hs+j,          ii) = haloRecvBufW(nUnpack*nGlob+iGlob);
      a(v,hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack*nGlob+iGlob);
    });
    nUnpack = nUnpack + n;
  }


  inline void haloUnpackN_y(Domain const &dom, real4d &a, int const n) {
    haloUnpackN_y_ext(dom, a, n, haloRecvBufS, haloRecvBufN, nUnpack);
  }
  inline void haloUnpackN_y_ext(Domain const &dom, real4d &a, int const n, real1d const &haloRecvBufS, real1d const &haloRecvBufN, int &nUnpack) {
    Kokkos::parallel_for( n*dom.nz*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, ii, i;
      unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
      int nGlob = dom.nz*hs*dom.nx;
      a(v,hs+k,          ii,hs+i) = haloRecvBufS(nUnpack*nGlob+iGlob);
      a(v,hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack*nGlob+iGlob);
    });
    nUnpack = nUnpack + n;
  }


  inline void haloExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    if (par.nproc_x > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( haloRecvBufW_cpu , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( haloRecvBufE_cpu , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

      #ifdef __NVCC__
        cudaMemcpyAsync( haloSendBufW_cpu , haloSendBufW.data() , nPack*dom.nz*dom.ny*hs*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaMemcpyAsync( haloSendBufE_cpu , haloSendBufE.data() , nPack*dom.nz*dom.ny*hs*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
      #endif

      //Send the data
      ierr = MPI_Isend( haloSendBufW_cpu , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( haloSendBufE_cpu , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      #ifdef __NVCC__
        cudaMemcpyAsync( haloRecvBufW.data() , haloRecvBufW_cpu , nPack*dom.nz*dom.ny*hs*sizeof(real) , cudaMemcpyHostToDevice );
        cudaMemcpyAsync( haloRecvBufE.data() , haloRecvBufE_cpu , nPack*dom.nz*dom.ny*hs*sizeof(real) , cudaMemcpyHostToDevice );
        cudaDeviceSynchronize();
      #endif

    } else {
      haloExchange_x_loc(dom, haloSendBufW, haloSendBufE, haloRecvBufW, haloRecvBufE);
    }
  }
  inline void haloExchange_x_loc(Domain const &dom, real1d &haloSendBufW, real1d &haloSendBufE, real1d &haloRecvBufW, real1d &haloRecvBufE) {
    Kokkos::parallel_for( nPack*dom.nz*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      haloRecvBufW(iGlob) = haloSendBufE(iGlob);
      haloRecvBufE(iGlob) = haloSendBufW(iGlob);
    });
  }


  inline void haloExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    if (par.nproc_y > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( haloRecvBufS_cpu , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( haloRecvBufN_cpu , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

      #ifdef __NVCC__
        cudaMemcpyAsync( haloSendBufS_cpu , haloSendBufS.data() , nPack*dom.nz*hs*dom.nx*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaMemcpyAsync( haloSendBufN_cpu , haloSendBufN.data() , nPack*dom.nz*hs*dom.nx*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
      #endif

      //Send the data
      ierr = MPI_Isend( haloSendBufS_cpu , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( haloSendBufN_cpu , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      #ifdef __NVCC__
        cudaMemcpyAsync( haloRecvBufS.data() , haloRecvBufS_cpu , nPack*dom.nz*hs*dom.nx*sizeof(real) , cudaMemcpyHostToDevice );
        cudaMemcpyAsync( haloRecvBufN.data() , haloRecvBufN_cpu , nPack*dom.nz*hs*dom.nx*sizeof(real) , cudaMemcpyHostToDevice );
        cudaDeviceSynchronize();
      #endif

    } else {
      haloExchange_y_loc(dom, haloSendBufS, haloSendBufN, haloRecvBufS, haloRecvBufN);
    }
  }
  inline void haloExchange_y_loc(Domain const &dom, real1d &haloSendBufS, real1d &haloSendBufN, real1d &haloRecvBufS, real1d &haloRecvBufN) {
    Kokkos::parallel_for( nPack*dom.nz*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      haloRecvBufS(iGlob) = haloSendBufN(iGlob);
      haloRecvBufN(iGlob) = haloSendBufS(iGlob);
    });
  }


  inline void edgePackN_x(Domain const &dom, real5d const &a, int const n) {
    edgePackN_x_ext(dom,a,n,edgeSendBufW,edgeSendBufE,nPack);
  }
  inline void edgePackN_x_ext(Domain const &dom, real5d const &a, int const n, real1d &edgeSendBufW, real1d &edgeSendBufE, int &nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j;
      unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
      int nGlob = dom.nz*dom.ny;
      edgeSendBufW(nPack*nGlob+iGlob) = a(v,1,k,j,0     );
      edgeSendBufE(nPack*nGlob+iGlob) = a(v,0,k,j,dom.nx);
    });
    nPack = nPack + n;
  }


  inline void edgePackN_y(Domain const &dom, real5d const &a, int const n) {
    edgePackN_y_ext(dom,a,n,edgeSendBufS,edgeSendBufN,nPack);
  }
  inline void edgePackN_y_ext(Domain const &dom, real5d const &a, int const n, real1d &edgeSendBufS, real1d &edgeSendBufN, int &nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, i;
      unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
      int nGlob = dom.nz*dom.nx;
      edgeSendBufS(nPack*nGlob+iGlob) = a(v,1,k,0     ,i);
      edgeSendBufN(nPack*nGlob+iGlob) = a(v,0,k,dom.ny,i);
    });
    nPack = nPack + n;
  }


  inline void edgeUnpackN_x(Domain const &dom, real5d &a, int const n) {
    edgeUnpackN_x_ext(dom, a, n, edgeRecvBufW, edgeRecvBufE, nUnpack);
  }
  inline void edgeUnpackN_x_ext(Domain const &dom, real5d &a, int const n, real1d const &edgeRecvBufW, real1d const &edgeRecvBufE, int &nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.nz*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, k, j;
      unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
      int nGlob = dom.nz*dom.ny;
      a(v,0,k,j,0     ) = edgeRecvBufW(nUnpack*nGlob+iGlob);
      a(v,1,k,j,dom.nx) = edgeRecvBufE(nUnpack*nGlob+iGlob);
    });
    nUnpack = nUnpack + n;
  }


  inline void edgeUnpackN_y(Domain const &dom, real5d &a, int const n) {
    edgeUnpackN_y_ext(dom, a, n, edgeRecvBufS, edgeRecvBufN, nUnpack);
  }
  inline void edgeUnpackN_y_ext(Domain const &dom, real5d &a, int const n, real1d const &edgeRecvBufS, real1d const &edgeRecvBufN, int &nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int k=0; k<dom.nz; k++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nz*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, k, i;
      unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
      int nGlob = dom.nz*dom.nx;
      a(v,0,k,0     ,i) = edgeRecvBufS(nUnpack*nGlob+iGlob);
      a(v,1,k,dom.ny,i) = edgeRecvBufN(nUnpack*nGlob+iGlob);
    });
    nUnpack = nUnpack + n;
  }


  inline void edgeExchange_x(Domain const &dom, Parallel const &par) {
    int ierr;

    if (par.nproc_x > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( edgeRecvBufW_cpu , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( edgeRecvBufE_cpu , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

      #ifdef __NVCC__
        cudaMemcpyAsync( edgeSendBufW_cpu , edgeSendBufW.data() , nPack*dom.nz*dom.ny*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaMemcpyAsync( edgeSendBufE_cpu , edgeSendBufE.data() , nPack*dom.nz*dom.ny*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
      #endif

      //Send the data
      ierr = MPI_Isend( edgeSendBufW_cpu , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( edgeSendBufE_cpu , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      #ifdef __NVCC__
        cudaMemcpyAsync( edgeRecvBufW.data() , edgeRecvBufW_cpu , nPack*dom.nz*dom.ny*sizeof(real) , cudaMemcpyHostToDevice );
        cudaMemcpyAsync( edgeRecvBufE.data() , edgeRecvBufE_cpu , nPack*dom.nz*dom.ny*sizeof(real) , cudaMemcpyHostToDevice );
        cudaDeviceSynchronize();
      #endif

    } else {
      edgeExchange_x_loc(dom, edgeSendBufW, edgeSendBufE, edgeRecvBufW, edgeRecvBufE);
    }
  }
  inline void edgeExchange_x_loc(Domain const &dom, real1d &edgeSendBufW, real1d &edgeSendBufE, real1d &edgeRecvBufW, real1d &edgeRecvBufE) {
    Kokkos::parallel_for( nPack*dom.nz*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      edgeRecvBufW(iGlob) = edgeSendBufE(iGlob);
      edgeRecvBufE(iGlob) = edgeSendBufW(iGlob);
    });
  }


  inline void edgeExchange_y(Domain const &dom, Parallel const &par) {
    int ierr;

    if (par.nproc_y > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( edgeRecvBufS_cpu , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( edgeRecvBufN_cpu , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

      #ifdef __NVCC__
        cudaMemcpyAsync( edgeSendBufS_cpu , edgeSendBufS.data() , nPack*dom.nz*dom.nx*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaMemcpyAsync( edgeSendBufN_cpu , edgeSendBufN.data() , nPack*dom.nz*dom.nx*sizeof(real) , cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
      #endif

      //Send the data
      ierr = MPI_Isend( edgeSendBufS_cpu , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( edgeSendBufN_cpu , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      #ifdef __NVCC__
        cudaMemcpyAsync( edgeRecvBufS.data() , edgeRecvBufS_cpu , nPack*dom.nz*dom.nx*sizeof(real) , cudaMemcpyHostToDevice );
        cudaMemcpyAsync( edgeRecvBufN.data() , edgeRecvBufN_cpu , nPack*dom.nz*dom.nx*sizeof(real) , cudaMemcpyHostToDevice );
        cudaDeviceSynchronize();
      #endif

    } else {
      edgeExchange_y_loc(dom, edgeSendBufS, edgeSendBufN, edgeRecvBufS, edgeRecvBufN);
    }
  }
  inline void edgeExchange_y_loc(Domain const &dom, real1d &edgeSendBufS, real1d &edgeSendBufN, real1d &edgeRecvBufS, real1d &edgeRecvBufN) {
    Kokkos::parallel_for( nPack*dom.nz*dom.nx, KOKKOS_LAMBDA (int iGlob) {
      edgeRecvBufS(iGlob) = edgeSendBufN(iGlob);
      edgeRecvBufN(iGlob) = edgeSendBufS(iGlob);
    });
  }

};

#endif
