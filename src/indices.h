// indices.h
//
// Header file for the indices module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _INDICES_H_
#define _INDICES_H_

inline int ij2f(int *ij, int fres[ND]) {
#if ND==2
  return ij[1]*fres[0] + ij[0];
#elif ND==3
  return ij[2]*fres[1]*fres[0] + ij[1]*fres[0] + ij[0];
#elif ND==4
  return ij[3]*fres[2]*fres[1]*fres[0] + ij[2]*fres[1]*fres[0] + ij[1]*fres[0] + ij[0];
#elif ND==6
  return ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] + ij[4]*fres[3]*fres[2]*fres[1]*fres[0] + ij[3]*fres[2]*fres[1]*fres[0] + ij[2]*fres[1]*fres[0] + ij[1]*fres[0] + ij[0];
#endif
}

inline void f2ij(int f, int *ij, int fres[ND]) {
  
#if ND==2
  ij[1] = f / fres[0];
  ij[0] = f - ij[1] * fres[0] ;
#elif ND==3
  ij[2] = f / (fres[1]*fres[0]);
  ij[1] = (f-ij[2]*fres[1]*fres[0]) / fres[0];
  ij[0] = f - ij[2]*fres[1]*fres[0] - ij[1] * fres[0];
#elif ND==4
  ij[3] = f / (fres[2]*fres[1]*fres[0]);
  ij[2] = (f-ij[3]*fres[2]*fres[1]*fres[0]) / (fres[1]*fres[0]);
  ij[1] = (f-ij[3]*fres[2]*fres[1]*fres[0]-ij[2]*fres[1]*fres[0]) / fres[0];
  ij[0] = f - ij[3]*fres[2]*fres[1]*fres[0] - ij[2]*fres[1]*fres[0] - ij[1] * fres[0];
#elif ND==6
  ij[5] = f / (fres[4]*fres[3]*fres[2]*fres[1]*fres[0]);
  ij[4] = (f-ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0]) / (fres[3]*fres[2]*fres[1]*fres[0]);
  ij[3] = (f-ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[4]*fres[3]*fres[2]*fres[1]*fres[0]) / (fres[2]*fres[1]*fres[0]);
  ij[2] = (f-ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[3]*fres[2]*fres[1]*fres[0]) / (fres[1]*fres[0]);
  ij[1] = (f-ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[3]*fres[2]*fres[1]*fres[0] - ij[2]*fres[1]*fres[0]) / fres[0];
  ij[0] = f-ij[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[4]*fres[3]*fres[2]*fres[1]*fres[0] - ij[3]*fres[2]*fres[1]*fres[0] - ij[2]*fres[1]*fres[0] - ij[1]*fres[0];
#endif
}

inline int ij2f(int *ij, int fres[ND], int d) {
#if ND==2
  return d*fres[1]*fres[0] + ij2f(ij,fres);
#elif ND==3
  return d*fres[2]*fres[1]*fres[0] + ij2f(ij,fres);
#elif ND==4
  return d*fres[3]*fres[2]*fres[1]*fres[0] + ij2f(ij,fres);
#elif ND==6
  return d*fres[5]*fres[4]*fres[3]*fres[2]*fres[1]*fres[0] + ij2f(ij,fres);
#endif
}





#endif


