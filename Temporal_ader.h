
#pragma once

#include "const.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

/* REQUIRED:
int  static constexpr nTimeDerivs = [# tendency time derivatives needed by the time stepping scheme];
bool static constexpr timeAvg     = [whether the spatial operator should return a time-averaged tendency];

static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

void init(std::string inFile):
    - Process input file with filename "inFile"
    - Allocate and initialize internal stuff

void timeStep( StateArr &state , real dt ): 
    - Perform a single time step

const char * getTemporalName():
    - Return the name and info about this temporal operator
*/
template <class Spatial> class Temporal_ader {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  typedef typename Spatial::StateTendArr  StateTendArr;
  typedef typename Spatial::StateArr      StateArr;
  typedef typename Spatial::TracerTendArr TracerTendArr;
  typedef typename Spatial::TracerArr     TracerArr;
  typedef typename Spatial::Location      Location;

  StateTendArr  stateTendArr;
  TracerTendArr tracerTendArr;

  Spatial spaceOp;
  
  void init(std::string inFile) {
    spaceOp.init(inFile);
    stateTendArr  = spaceOp.createStateTendArr ();
    tracerTendArr = spaceOp.createTracerTendArr();
  }


  void timeStep( StateArr &stateArr , TracerArr &tracerArr , real dt ) {
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      int constexpr hs       = Spatial::hs;
      int constexpr numState = Spatial::numState;
      int numTracers = spaceOp.numTracers;
      int nz = spaceOp.nz;
      int ny = spaceOp.ny;
      int nx = spaceOp.nx;
      auto &stateTendArr  = this->stateTendArr ;
      auto &tracerTendArr = this->tracerTendArr;

      real dtloc = dt;
      spaceOp.computeTendencies( stateArr , stateTendArr , tracerArr , tracerTendArr , dtloc , spl );

      parallel_for( SimpleBounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
        stateArr(l,hs+k,hs+j,hs+i) += dtloc * stateTendArr(l,k,j,i);
      });

      parallel_for( SimpleBounds<4>(numTracers,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
        tracerArr(l,hs+k,hs+j,hs+i) += dtloc * tracerTendArr(l,k,j,i);
      });

    }
  }


  void finalize(StateArr &state, TracerArr &tracers) {
    spaceOp.finalize( state , tracers );
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

