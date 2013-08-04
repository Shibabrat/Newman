// boundary.h
//
// Header file for the boundary module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "velocity.h"
#include "macros.h"
#include "errors.h"
#include "mapcoord.h"
#include "indices.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#define MAXNUMISLANDS 100
#define MAXNUMSEGMENTS 5000

extern int  myrank;

extern double FBP_mu;
extern double FBP_a;
extern double FBP_omega;
extern double FBP_theta0;
extern double FBP_mu3;
extern double FBP_sunbdy2;
extern double FBP_moonbdy2;
extern double FBP_earthbdy2;

extern int TBP_Ordinate;

extern double **DataNest_Min;
extern double **DataNest_Max;
extern vector<int> Nest_List;

extern int  Velocity_Format;
extern int Data_Periodic[ND];
extern double Data_Period[ND];
extern double Data_Min[ND];
extern double Data_Max[ND];
extern double Data_Delta[ND];
extern int Data_Res[ND];
extern int    Data_HeaderSize;
extern int  MapCoord;

extern int Boundary_Method;
extern double Boundary_MaskValue;
extern int Boundary_NumLists;
extern char Boundary_Input[SHORTSTRING];

extern char Path_Output[LONGSTRING];

struct Boundary_List {
  int    numsegments;
  double **X;
  double *limits;
};

// Function Prototypes:
void CleanUpBoundary();
void ReadInBoundaryData(int format);
void ReadInTecBoundaryData();
void GetBoundaryMask();
void WriteOutTecBoundary();
int TestOutsideAnalyticalBoundary(double *point, double t);
int TestOutsideStaircase(double *point);
int TestOutsidePolygon(double *point); 
int LongTestOutsidePolygon(int jbdy,double *point);
int QuickTestOutsidePolygon(int jbdy,double *point);
void Setup_Boundary();

inline int TestOutsideDomain(double *point, double t) {
  
  /*  
  return 0 -> still in flow domain
  return 1 -> Still in data rectangle but outside flow. 
  Must decide whether to enforce no flow. (Not Yet Supported) 
  return 2 -> left data rectangle entirely.
  */
  
  
  if(Velocity_Format == 2) {
    for(int d=0;d<ND;++d) {
      if(!Data_Periodic[d]){
        if( point[d] < Data_Min[d] || point[d] > Data_Max[d]) {
          return(1);
        }
      }
    }
  }

  if(Velocity_Format == 3) {
    for(int d=0;d<ND;++d) {
        if( point[d] < Data_Min[d] || point[d] > Data_Max[d]) {
          return(1);
      }
    }
  }

  
  if(Boundary_Method) {
    switch ( Boundary_Method ) {
      
      case 1 : //Use staircase boundary
        if(TestOutsideStaircase(point))
          return(1);
        break;
        
#if ND==2  
      case 2 : //Use polygon boundary
        if(TestOutsidePolygon(point))
          return(1);
        break;
#endif
      case 3 : //Use analytical boundary
        if(TestOutsideAnalyticalBoundary(point,t))
          return(1);
      break;

    }
  }
    
  return(0);
}





#endif


