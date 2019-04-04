
#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "const.h"

class Parallel {

public:

  ulong nTasks;
  ulong myTask;
  ulong px;
  ulong py;
  ulong ix;
  ulong iy;
  int masterTask;
};

#endif
