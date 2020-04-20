
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "mpi.h"
#include "Domain.h"
#include "Parallel.h"


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

  realHost1d haloSendBufS_host;
  realHost1d haloSendBufN_host;
  realHost1d haloSendBufW_host;
  realHost1d haloSendBufE_host;
  realHost1d haloRecvBufS_host;
  realHost1d haloRecvBufN_host;
  realHost1d haloRecvBufW_host;
  realHost1d haloRecvBufE_host;

  realHost1d edgeRecvBufE_host;
  realHost1d edgeRecvBufW_host;
  realHost1d edgeSendBufE_host;
  realHost1d edgeSendBufW_host;
  realHost1d edgeRecvBufN_host;
  realHost1d edgeRecvBufS_host;
  realHost1d edgeSendBufN_host;
  realHost1d edgeSendBufS_host;

public:


  void allocate(Domain &dom);


  void haloInit();


  void haloPackN_x(Domain const &dom, real4d const &a, int const n);


  void haloPackN_y(Domain const &dom, real4d const &a, int const n);


  void haloUnpackN_x(Domain const &dom, real4d &a, int const n);


  void haloUnpackN_y(Domain const &dom, real4d &a, int const n);


  void haloExchange_x(Domain const &dom, Parallel const &par);
  void haloExchange_x_loc(Domain const &dom, real1d &haloSendBufW, real1d &haloSendBufE, real1d &haloRecvBufW, real1d &haloRecvBufE);


  void haloExchange_y(Domain const &dom, Parallel const &par);
  void haloExchange_y_loc(Domain const &dom, real1d &haloSendBufS, real1d &haloSendBufN, real1d &haloRecvBufS, real1d &haloRecvBufN);


  void edgePackN_x(Domain const &dom, real5d const &a, int const n);


  void edgePackN_y(Domain const &dom, real5d const &a, int const n);


  void edgeUnpackN_x(Domain const &dom, real5d &a, int const n);


  void edgeUnpackN_y(Domain const &dom, real5d &a, int const n);


  void edgeExchange_x(Domain const &dom, Parallel const &par);
  void edgeExchange_x_loc(Domain const &dom, real1d &edgeSendBufW, real1d &edgeSendBufE, real1d &edgeRecvBufW, real1d &edgeRecvBufE);


  void edgeExchange_y(Domain const &dom, Parallel const &par);
  void edgeExchange_y_loc(Domain const &dom, real1d &edgeSendBufS, real1d &edgeSendBufN, real1d &edgeRecvBufS, real1d &edgeRecvBufN);

};

#endif
