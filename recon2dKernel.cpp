
#include "YAKL.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"

#ifndef WENO
  #define WENO true;
#endif

int  constexpr nx       = 200;
int  constexpr ny       = 200;
int  constexpr nz       = 100;
bool constexpr doWeno   = WENO;



int main() {
  yakl::init();

  std::cout << "ord:  " << ord    << "\n";
  std::cout << "tord: " << tord   << "\n";
  std::cout << "weno: " << doWeno << "\n";

  TransformMatrices<real> trans;

  SArray<real,2,ord,tord> to_gll;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> wenoIdl;
  real wenoSigma;

  if (doWeno) {
    trans.coefs_to_gll_lower( to_gll );
  } else {
    trans.sten_to_gll_lower( to_gll );
  }
  trans.weno_sten_to_coefs( wenoRecon );
  wenoSetIdealSigma(wenoIdl,wenoSigma);
  
  real4d state("state",numState,nz+2*hs,ny+2*hs,nx+2*hs);
  real6d stateGLL("stateGLL",numState,nz,ny,nx,tord,tord);

  parallel_for( Bounds<4>(numState,nz+2*hs,ny+2*hs,nx+2*hs) , YAKL_LAMBDA (int v, int k, int j, int i) {
    state(v,k,j,i) = 1;
  });

  if (doWeno) {

    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      for (int l=0; l<numState; l++) {
        SArray<real,1,ord> avg;
        SArray<real,1,hs+2> wenoWts;
        SArray<real,1,ord> coefs;
        SArray<real,2,ord,tord> glltmp2d;
        //////////////////////////////////////////////////////////////////////////////////////
        // X-direction
        //////////////////////////////////////////////////////////////////////////////////////
        // Compute y-direction average
        for (int ii=0; ii<ord; ii++) { avg(ii) = 0; }
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            avg(ii) += state(l,hs+k,j+jj,i+ii);
          }
        }
        for (int ii=0; ii<ord; ii++) { avg(ii) /= ord; }
        // Compute weno weights on y-direction average
        weno_recon_and_weights( wenoRecon , avg , wenoIdl , wenoSigma , wenoWts );
        // Apply WENO weights to each jj
        for (int jj=0; jj<ord; jj++) {
          SArray<real,1,ord> sten1d;
          for (int ii=0; ii<ord; ii++) { sten1d(ii) = state(l,hs+k,j+jj,i+ii); }
          weno_recon_and_apply( wenoRecon , sten1d , wenoIdl , wenoWts , coefs );
          for (int ii=0; ii < tord; ii++) {
            real tmp = 0;
            for (int s=0; s < ord; s++) {
              tmp += to_gll(s,ii) * coefs(s);
            }
            glltmp2d(jj,ii) = tmp;
          }
        }

        //////////////////////////////////////////////////////////////////////////////////////
        // Y-direction
        //////////////////////////////////////////////////////////////////////////////////////
        // Compute x-direction average
        for (int jj=0; jj<ord; jj++) { avg(jj) = 0; }
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<tord; ii++) {
            avg(jj) += glltmp2d(jj,ii);
          }
        }
        for (int jj=0; jj<ord; jj++) { avg(jj) /= tord; }
        // Compute weno weights on x-direction average
        weno_recon_and_weights( wenoRecon , avg , wenoIdl , wenoSigma , wenoWts );
        // Apply WENO weights to each ii
        for (int ii=0; ii<tord; ii++) {
          SArray<real,1,ord> sten1d;
          for (int jj=0; jj<ord; jj++) { sten1d(jj) = glltmp2d(jj,ii); }
          weno_recon_and_apply( wenoRecon , sten1d , wenoIdl , wenoWts , coefs );
          for (int jj=0; jj < tord; jj++) {
            real tmp = 0;
            for (int s=0; s < ord; s++) {
              tmp += to_gll(s,jj) * coefs(s);
            }
            stateGLL(l,k,j,i,jj,ii) = tmp;
          }
        }
      }
    });

  } else {

    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      for (int l=0; l<numState; l++) {
        SArray<real,2,ord,tord> glltmp2d;
        //////////////////////////////////////////////////////////////////////////////////////
        // X-direction
        //////////////////////////////////////////////////////////////////////////////////////
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii < tord; ii++) {
            real tmp = 0;
            for (int s=0; s < ord; s++) {
              tmp += to_gll(s,ii) * state(l,hs+k,j+jj,i+s);
            }
            glltmp2d(jj,ii) = tmp;
          }
        }

        //////////////////////////////////////////////////////////////////////////////////////
        // Y-direction
        //////////////////////////////////////////////////////////////////////////////////////
        for (int ii=0; ii<tord; ii++) {
          for (int jj=0; jj < tord; jj++) {
            real tmp = 0;
            for (int s=0; s < ord; s++) {
              tmp += to_gll(s,jj) * glltmp2d(s,ii);
            }
            stateGLL(l,k,j,i,jj,ii) = tmp;
          }
        }
      }
    });

  }

  yakl::finalize();
}



