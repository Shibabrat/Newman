// mapcoord.cpp
//
// mapcoord module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "mapcoord.h"

void MapCoordToIntegrationSpace(double *X, double t) {

#if ND==2
  /*double x = X[0];
  double y = X[1];
  
  X[0] = hypot(x,y);
  X[1] = atan2(y,x);*/
  
  X[0] = X[0]*sqrt(DNA_N);
  X[1] = X[1]*sqrt(DNA_N);

#elif ND==4
  if (TBP_Ordinate == 0) {
  /* Map ( x , xdot , E , y )  to  ( x , y , xdot , ydot ) */
  double tmp = X[3];
  X[3] = TBP_VDirection*sqrt(2.0*(X[2]-0.5*X[1]*X[1]+(0.5*X[0]*X[0]+0.5*X[3]*X[3]+0.5*TBP_mu-0.5*TBP_mu*TBP_mu+TBP_mu/sqrt(X[3]*X[3]+(X[0]-1.0+TBP_mu)*(X[0]-1.0+TBP_mu))+(1.0-TBP_mu)/sqrt(X[3]*X[3]+(X[0]+TBP_mu)*(X[0]+TBP_mu)))/(1.0+TBP_ecc*cos(t))));
  X[2] = X[1];
  X[1] = tmp;
  }

  if (TBP_Ordinate == 1) {
  /* Map ( y , ydot , E , x )  to ( x , y , xdot , ydot ) */
  double tmp = X[3];
  X[2] = TBP_VDirection*sqrt(2.0*(X[2]-0.5*X[1]*X[1]+(0.5*X[3]*X[3]+0.5*X[0]*X[0]+0.5*TBP_mu-0.5*TBP_mu*TBP_mu+TBP_mu/sqrt(X[0]*X[0]+(X[3]-1.0+TBP_mu)*
(X[3]-1.0+TBP_mu))+(1.0-TBP_mu)/sqrt(X[0]*X[0]+(X[3]+TBP_mu)*(X[3]+TBP_mu)))/(1.0+TBP_ecc*cos(t))));
  X[3] = X[1];
  X[1] = X[0];
  X[0] = tmp;  
  }

#elif ND==6
  if (TBP_Ordinate == 0) {
  /* Map ( x , xdot , E , y , z, zdot )  to  ( x , y , z, xdot , ydot, zdot ) */
  double x = X[0];
  double vx = X[1];
  double E = X[2];
  double y = X[3];
  double z = X[4];
  double vz = X[5];
  double den = 1.0/(1.0+TBP_ecc*cos(t));
  double r1den = 1.0/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z);
  double r2den = 1.0/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z);
  double vy = TBP_VDirection*sqrt(2.0*(    E - 0.5*(vx*vx+vz*vz) + den*(0.5*(x*x+y*y)+0.5*TBP_mu*(1.0-TBP_mu)+TBP_mu*r2den+(1.0-TBP_mu)*r1den-0.5*TBP_ecc*z*z*cos(t))      ));
  X[0] = x;
  X[1] = y;
  X[2] = z;
  X[3] = vx;
  X[4] = vy;
  X[5] = vz;
  }

  if (TBP_Ordinate == 1) {
  /* Map ( y , ydot , E , x , z, zdot )  to  ( x , y , z, xdot , ydot, zdot ) */
  double y = X[0];
  double vy = X[1];
  double E = X[2];
  double x = X[3];
  double z = X[4];
  double vz = X[5];
  double den = 1.0/(1.0+TBP_ecc*cos(t));
  double r1den = 1.0/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z);
  double r2den = 1.0/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z);
  double vx = TBP_VDirection*sqrt(2.0*(    E - 0.5*(vy*vy+vz*vz) + den*(0.5*(x*x+y*y)+0.5*TBP_mu*(1.0-TBP_mu)+TBP_mu*r2den+(1.0-TBP_mu)*r1den-0.5*TBP_ecc*z*z*cos(t))      ));
  X[0] = x;
  X[1] = y;
  X[2] = z;
  X[3] = vx;
  X[4] = vy;
  X[5] = vz;
  }

#endif
  
}

void MapIntegrationSpaceToCoord(double *X, double t) {
  
#if ND==2
  /*double r = X[0];
  double theta = X[1];
  
  X[0] = r*cos(theta);
  X[1] = r*sin(theta);*/
  
  X[0] = X[0]/sqrt(DNA_N);
  X[1] = X[1]/sqrt(DNA_N);

#elif ND==4
  if (TBP_Ordinate == 0) {
    /* Map ( x , y , xdot , ydot ) to ( x , xdot , E , y ) */
    double tmp = X[2];
    X[2]=0.5*(X[2]*X[2]+X[3]*X[3])+(-0.5*(X[0]*X[0]+X[1]*X[1])-(1-TBP_mu)/sqrt((X[0]+TBP_mu)*(X[0]+TBP_mu)+X[1]*X[1])-TBP_mu/sqrt((X[0]-1.0+TBP_mu)*(X[0]-1.0+TBP_mu)+X[1]*X[1])-0.5*TBP_mu*(1.0-TBP_mu))/(1.0+TBP_ecc*cos(t));
    X[3]=X[1];
    X[1]=tmp;
  }

  if (TBP_Ordinate == 1) {
  /* Map ( x , y , xdot , ydot ) to ( y , ydot , E , x ) */
  double tmp = X[3];
  X[2]=0.5*(X[2]*X[2]+X[3]*X[3])+(-0.5*(X[0]*X[0]+X[1]*X[1])-(1-TBP_mu)/sqrt((X[0]+TBP_mu)*(X[0]+TBP_mu)+X[1]*X[1])-TBP_mu/sqrt((X[0]-1.0+TBP_mu)*(X[0]-1.0+TBP_mu)+X[1]*X[1])-0.5*TBP_mu*(1.0-TBP_mu))/(1.0+TBP_ecc*cos(t));
  X[3]=X[0];
  X[0]=X[1];
  X[1]=tmp;
  }

#elif ND==6
  if (TBP_Ordinate == 0) {
  /* Map (x, y, z, xdot, ydot, zdot) to ( x, xdot, E, y, z, zdot ) */
  double x = X[0];
  double y = X[1];
  double z = X[2];
  double vx = X[3];
  double vy = X[4];
  double vz = X[5];
  double E = 0.5*(vx*vx+vy*vy+vz*vz)+(-0.5*(x*x+y*y-z*z*TBP_ecc*cos(t))-(1.0-TBP_mu)/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z)-TBP_mu/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z)-0.5*TBP_mu*(1.0-TBP_mu))/(1.0+TBP_ecc*cos(t));
  X[0] = x;
  X[1] = vx;
  X[2] = E;
  X[3] = y;
  X[4] = z;
  X[5] = vz;
  }

  if (TBP_Ordinate == 1) {
  /* Map (x, y, z, xdot, ydot, zdot)  to ( y , ydot , E , x, z, zdot ) */
  double x = X[0];
  double y = X[1];
  double z = X[2];
  double vx = X[3];
  double vy = X[4];
  double vz = X[5];
  double E = 0.5*(vx*vx+vy*vy+vz*vz)+(-0.5*(x*x+y*y-z*z*TBP_ecc*cos(t))-(1.0-TBP_mu)/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z)-TBP_mu/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z)-0.5*TBP_mu*(1.0-TBP_mu))/(1.0+TBP_ecc*cos(t));
  X[0] = y;
  X[1] = vy;
  X[2] = E;
  X[3] = x;
  X[4] = z;
  X[5] = vz;
  }

#endif
  
  // Planar: E=1/2*(vx^2+vy^2)+(-1/2*(x^2+y^2)-(1-mu)/sqrt((x+mu)^2+y^2)-mu/sqrt((x-1+mu)^2+y^2)-1/2*mu*(1-mu))/(1+ecc*cos(t));
  // Spatial: E=1/2*(vx^2+vy^2+vz^2)+(-1/2*(x^2+y^2-z^2*ecc*cos(t))-(1-mu)/sqrt((x+mu)^2+y^2+z^2)-mu/sqrt((x-1+mu)^2+y^2+z^2)-1/2*mu*(1-mu))/(1+ecc*cos(t));
  
}





