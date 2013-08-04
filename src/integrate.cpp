// integrate.cpp
//
// integrate module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include <cmath>
#include "integrate.h"
#include "errors.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>



double params;


const gsl_odeiv_step_type *T;
gsl_odeiv_step *s;
gsl_odeiv_control *con;
gsl_odeiv_evolve *e;

int func (double t, const double X[], double dXdt[], void *params) {
  
  GetVelocity(t, (double *)X, dXdt);
  
  return GSL_SUCCESS;
}

int jacobian (double t, const double X[], double *dfdX, double dfdt[], void *params)
{
  
  GetCodedJacobian( t , (double *)X , dfdX , (double *)dfdt );
  return GSL_SUCCESS;
  
}

gsl_odeiv_system sys = {func, jacobian, ND , &params};

void SetUpIntegrator() {
  
  switch (Int_Method) {
    case 0:
      T = gsl_odeiv_step_rkf45;
      s = gsl_odeiv_step_alloc(T, ND);
      if(myrank == 0) cout << "Intializing GSL Adaptive step size RK45 integrator: " << gsl_odeiv_step_name(s) << endl;
        con = gsl_odeiv_control_y_new (Int_AbsTol,Int_RelTol);
      e = gsl_odeiv_evolve_alloc (ND);
      break;
    case 1:
      T = gsl_odeiv_step_rk4;
      s = gsl_odeiv_step_alloc(T, ND);
      if(myrank == 0) cout << "Intializing GSL Fixed step size RK 4 integrator: " << gsl_odeiv_step_name(s) << endl;
        break;
    case 2:
      if(myrank == 0) cout << "Using Euler integrator:" << endl;
      break;
    case 3:
      if(myrank == 0) cout << "Using Adaptive step size RK45 integrator" << endl;
      break;
  }
  return;
  
}

void CleanUpIntegrator() {
  
  switch (Int_Method) {
    case 0:
      gsl_odeiv_evolve_free(e);
      gsl_odeiv_control_free(con);
      gsl_odeiv_step_free(s);
      break;
    case 1:
      gsl_odeiv_step_free(s);
      break;
    case 2:
      break;
      
  }
  return;
  
}





int IntGSLfixed(double *X, double tstart, double tend) {
  
  
  double currt=tstart;
  double h = SIGN(Int_TimeStep, tend - tstart);
  
  double dydt_in[ND], dydt_out[ND];
  GetVelocity(currt,X,dydt_in);  
  double xerr[ND];
  while( (currt - tend) * (tend - tstart) < 0.0 ) {
    
    if ((currt + h - tend) * (currt + h - tstart) > 0.0) {
      h = tend - currt;  
    }
    
    int status = gsl_odeiv_step_apply (s, currt, h, X, xerr, dydt_in, dydt_out, &sys);
    if(status != GSL_SUCCESS)
      FatalError(gsl_strerror(status));
    
    for(int d=0;d<ND;++d) {
      dydt_in[d] = dydt_out[d];
      if(!isfinite(X[d])) {
        cout << "Integrator blow up." << endl;
        exit(1);
      }
    }
    currt += h;
    
  }
  
  
  return 0;
}


int IntGSLadapt(double *X, double tstart, double tend) {
  
  double currt=tstart;
  double h = SIGN(Int_TimeStep, tend - tstart);
  
  while( (currt - tend) * (tend - tstart) < 0.0 ) {
    
    int status = gsl_odeiv_evolve_apply (e, con, s, &sys, &currt, tend, &h, X);
    
    if(status != GSL_SUCCESS)
      FatalError(gsl_strerror(status));
  }
  
  return 0;
}

void IntEuler(double *X, double tstart, double tend) {
  /***** Advect points using simple Euler *****/
  
  double dXdt[ND];  
  double h = SIGN(Int_TimeStep, tend - tstart);
  double currt = tstart;
  
  while((currt - tend) * (tend - tstart) < 0.0 ) {
    if ((currt + h - tend) * (currt + h - tstart) > 0.0) {
      h = tend - currt;  
    }
    
    GetVelocity(currt, X, dXdt); 
    
    for(int d=0;d<ND;++d)
      X[d] += dXdt[d] * h;
    
    currt += h;
  }
}


void IntAdaptRK45(double *X, double tstart, double tend) {
  /***** Advect points using RKF45 and bicubic spatial interpolation *****
  *     temporal interpolation for RKF45 midpoints and self-adjusting   *
  *     step size is Nth (nd) order Lagrange polynomial                 */
  
  double k1[ND],k2[ND],k3[ND],k4[ND],k5[ND],k6[ND];
  double a3=3./32., b3=9./32., a4=1932./2197., b4=-7200./2197., c4=7296./2197.,
    a5=439./216., b5=-8., c5=3680./513., d5=-845./4104.,
    a6=-8./27., b6=2., c6=-3544./2565., d6=1859./4104., e6=-11./40.,
    a7=1./360., b7=-128./4275., c7=-2197./75240., d7=1./50., e7=2./55.,
    a8=25./216., b8=1408./2565., c8=2197./4104., d8=-1./5.;
  
  /* b7 is different from MANGEN */
  double errmax=0.0;
  double currtp[ND], err[ND], ss, tcurrtp;
  
  double h = SIGN(Int_TimeStep, tend - tstart);
  
  double currt = tstart;
  while((currt - tend) * (tend - tstart) < 0.0 ) {
    if ((currt + h - tend) * (currt + h - tstart) > 0.0) {
      h = tend - currt;
    }
    
    /* RK step 1 */
    GetVelocity(currt, X, k1);
    
    for(int d=0;d<ND;++d)
      currtp[d] = X[d] + 0.25 * k1[d] * h;
    
    /* RK step 2 */
    tcurrtp = currt + 0.25 * h;
    GetVelocity(tcurrtp, currtp, k2);
    for(int d=0;d<ND;++d)
      currtp[d] = X[d] + (a3 * k1[d] + b3 * k2[d]) * h;
    
    /* RK step 3 */
    tcurrtp = currt + 0.375 * h;
    GetVelocity(tcurrtp, currtp, k3);
    for(int d=0;d<ND;++d)
      currtp[d] = X[d] + (a4 * k1[d] + b4 * k2[d] + c4 * k3[d]) * h;
    
    /* RK step 4 */
    tcurrtp = currt + 12.0 * h / 13.0;
    GetVelocity(tcurrtp, currtp, k4);
    for(int d=0;d<ND;++d)
      currtp[d] = X[d] + (a5 * k1[d] + b5 * k2[d] + c5 * k3[d] + d5 * k4[d]) * h;
    
    /* RK step 5 */
    tcurrtp = currt + h;
    GetVelocity(tcurrtp, currtp, k5);
    for(int d=0;d<ND;++d)
      currtp[d] = X[d] + (a6 * k1[d] + b6 * k2[d] + c6 * k3[d] + d6 * k4[d] + e6 * k5[d]) * h;
    
    tcurrtp = currt + 0.5 * h;
    GetVelocity(tcurrtp, currtp, k6);
    for(int d=0;d<ND;++d)
      err[d] = fabs(a7 * k1[d] + b7 * k3[d] + c7 * k4[d] + d7 * k5[d] + e7 * k6[d]);
    
    errmax = 0.0;  /* This is a very important line */
    for(int d=0;d<ND;++d)
      errmax = fmax(errmax,err[d]);
    
    if (errmax == 0.)
      ss = 0.;
    else
      ss = pow(fabs(Int_AbsTol * h / (2.0 * errmax)), 0.25);
      
    if (errmax < Int_AbsTol ) {  /* accept approximation */
      for(int d=0;d<ND;++d)
        X[d] += (a8 * k1[d] + b8 * k3[d] + c8 * k4[d] + d8 * k5[d]) * h;
      
      currt = currt + h;
     
      /*** The below formulas for changing step size are not standard, google Runge-Kutta-Felhburg ***/
      if (ss > 1.5) { /* increase step size */
        h = SIGN((fabs(h * 2.0) > Int_MaxTimeStep)?Int_MaxTimeStep:h*2.0, tend - tstart);
      }
    }
    else if (fabs(h) <= Int_MinTimeStep) { /* accept approximation */
      for(int d=0;d<ND;++d)
        X[d] += (a8 * k1[d] + b8 * k3[d] + c8 * k4[d] + d8 * k5[d]) * h;

      currt += h;
    }
    else { /* do not accept */
      h = SIGN(((fabs(h / 2.0) < Int_MinTimeStep)?Int_MinTimeStep:h/2.0) , tend - tstart);
    }
  }
}





