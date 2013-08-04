// mapcoord.h
//
// Header file for the coordmap module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _MAPCOORD_H_
#define _MAPCOORD_H_

#include <cmath>
#include <iostream>
using namespace std;
#include "macros.h"

extern double TBP_mu;
extern double TBP_ecc;
extern int    TBP_Ordinate;
extern int    TBP_VDirection;
extern int DNA_N;

// Function Prototypes:
void MapCoordToIntegrationSpace(double *X, double t);
void MapIntegrationSpaceToCoord(double *X, double t);

#endif


