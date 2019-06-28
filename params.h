
#ifndef _PARAMS_H_
#define _PARAMS_H_

int const TIME_ADER   = 1;
int const TIME_SSPRK3 = 2;

int const EQN_THETA_CONS = 1;
int const EQN_THETA_PRIM = 2;

int const DATA_INIT_COLLISION = 1;
int const DATA_INIT_THERMAL   = 2;
int const DATA_INIT_STRAKA    = 3;

real outFreq;
int  timeMethod;
int  strakaVisc;

#endif
