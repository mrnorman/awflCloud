
#include "YAKL.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"

int  constexpr nx       = 50;
int  constexpr ny       = 50;
int  constexpr nz       = 25;
bool constexpr doWeno   = true;



int main() {
  yakl::init();

  TransformMatrices<real> trans;
  real dx = 20000/nx;
  real dy = 20000/ny;
  real dz = 20000/nz;
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

  parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
    for (int l=0; l<numState; l++) {
      if (doWeno) {
        SArray<real,2,ord,ord> stencil;
        SArray<real,1,ord> avg;
        SArray<real,1,hs+2> wenoWts;
        SArray<real,1,ord> coefs;
        SArray<real,1,tord> gll;
        SArray<real,2,ord,tord> glltmp2d;
        //////////////////////////////////////////////////////////////////////////////////////
        // X-direction
        //////////////////////////////////////////////////////////////////////////////////////
        // Load stencil
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            stencil(jj,ii) = state(l,hs+k,j+jj,i+ii);
          }
        }
        // Compute y-direction average
        for (int ii=0; ii<ord; ii++) { avg(ii) = 0; }
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            avg(ii) += stencil(jj,ii);
          }
        }
        for (int ii=0; ii<ord; ii++) { avg(ii) /= ord; }
        // Compute weno weights on y-direction average
        weno_recon_and_weights( wenoRecon , avg , wenoIdl , wenoSigma , wenoWts );
        // Apply WENO weights to each jj
        for (int jj=0; jj<ord; jj++) {
          SArray<real,1,ord> sten1d;
          for (int ii=0; ii<ord; ii++) { sten1d(ii) /= stencil(jj,ii); }
          weno_recon_and_apply( wenoRecon , sten1d , wenoIdl , wenoWts , coefs );
          gll = to_gll * coefs;
          for (int ii=0; ii<tord; ii++) { glltmp2d(jj,ii) = gll(ii); }
        }

        //////////////////////////////////////////////////////////////////////////////////////
        // Y-direction
        //////////////////////////////////////////////////////////////////////////////////////
        // Compute x-direction average
        for (int jj=0; jj<ord; jj++) { avg(jj) = 0; }
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<tord; ii++) {
            avg(jj) += stencil(jj,ii);
          }
        }
        for (int jj=0; jj<ord; jj++) { avg(jj) /= tord; }
        // Compute weno weights on y-direction average
        weno_recon_and_weights( wenoRecon , avg , wenoIdl , wenoSigma , wenoWts );
        // Apply WENO weights to each jj
        for (int ii=0; ii<tord; ii++) {
          SArray<real,1,ord> sten1d;
          for (int jj=0; jj<ord; jj++) { sten1d(jj) /= stencil(jj,ii); }
          weno_recon_and_apply( wenoRecon , sten1d , wenoIdl , wenoWts , coefs );
          gll = to_gll * coefs;
          for (int jj=0; jj<tord; jj++) { stateGLL(l,k,j,i,jj,ii) = gll(jj); }
        }
      } else {
      }
    }
  });

  yakl::finalize();
}



