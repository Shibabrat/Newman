// integrate.h
//
// Header file for integrate module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#ifndef _INTEGRATE_H_
#define _INTEGRATE_H_

#include "macros.h"
#include "velocity.h"

extern int    Int_Method;
extern double Int_TimeStep;
extern double Int_MinTimeStep;
extern double Int_MaxTimeStep;
extern double Int_AbsTol;
extern double Int_RelTol;
extern int    Int_UseJacobian;

void SetUpIntegrator();
void CleanUpIntegrator();
int IntGSLadapt(double *X, double tstart, double tend);
int IntGSLfixed(double *X, double tstart, double tend);
void IntAdaptRK45(double *X, double tstart, double tend);
void IntEuler(double *X, double tstart, double tend);
int func (double t, double y[], double f[], void *params);
int jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params);

inline void Integrate(double *X, double tstart, double tend) {

  if(MapCoord) MapCoordToIntegrationSpace(X,tstart);

  switch (Int_Method) {
    case 0:
      IntGSLadapt(X, tstart, tend);
      break;
    case 1:
      IntGSLfixed(X, tstart, tend);
      break;
    case 2:
      IntEuler(X, tstart, tend);
      break;
    case 3:
      IntAdaptRK45(X, tstart, tend);
      break;
    default:
      IntGSLadapt(X, tstart, tend);
  }
  
  if(MapCoord) MapIntegrationSpaceToCoord(X,tend);
  
  return;
}

#endif


