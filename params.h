
#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "YAKL.h"

int const TIME_ADER   = 1;
int const TIME_SSPRK3 = 2;

real outFreq;
int  doWeno;
int  timeMethod;

yakl::Launcher launcher(yakl::targetCUDA);

#endif
