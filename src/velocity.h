// velocity.cpp
//
// velocity module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#ifndef _VELOCITY_H_
#define _VELOCITY_H_

#include <vector>
#include "boundary.h"
#include "eqnparser.h"
#include "indices.h"

extern vector<int> Time_Origin;

extern char Plot_Work[LONGSTRING];
extern char Plot_Scratch[LONGSTRING];
extern char Plot_Output[LONGSTRING];

extern vector<double>  DNA_q0;
extern double **Utransform;

extern int DNA_N;
extern double ABC_Amplitude;

extern double Velocity_Null[ND];

extern  int    Query_Velocity;
extern  double Query_X[ND];

extern vector<int> Nest_List;

extern int    Velocity_Format;

extern int    Cell_NumVerts;
extern double *Cell_Areas;
extern int    Cell_Mask[ND];

extern int    Time_Direction;
extern double datatime1;
extern int    Data_TPeriodic;
extern double Data_TPeriod;
extern int    Atmos_Set;
extern double Atmos_Radius;
extern double Atmos_AltScale;
extern int    Velocity_Format;
extern double **Data_Array;
extern int    Data_BlockSize;

extern double Data_Max[ND];
extern int    Data_Periodic[ND];
extern double Data_Period[ND];
extern double Data_TDelta;
extern double Data_Min[ND];
extern double Data_Delta[ND];
extern int    Data_Res[ND];
extern int    Data_NonUniformGrid[ND];
extern double *Data_Grid[ND];
extern int *Boundary_Mask;

extern double ***DataNest_Array;
extern double **DataNest_Min;
extern double **DataNest_Max;
extern double **DataNest_Delta;
extern int **DataNest_Res;

extern double Plot_Min[ND];
extern double Plot_Max[ND];
extern int    Plot_Res[ND];
extern double Plot_Delta[ND];
extern int    Plot_BlockSize;

extern char Path_Work[LONGSTRING];
extern char   Plot_OutFile[SHORTSTRING];
extern char   Plot_Output[LONGSTRING];

extern char AnalyticEq[ND][LONGSTRING];

//Function prototypes
void  Hat2Reg(vector<double>& qhat, vector<double>& qreg);
void  Reg2Hat(vector<double>& qreg, vector<double>& qhat);
void  Lift(double t, double *q, vector<double>& qfull);
double  LinSol(double t,int w);
void GenerateU();
void GetVelocity(double t,double *X,double *dXdt);
void GetCartesianVelocity(double t, double *X, double *dXdt);
void GetCodedVelocity(double t,double *X,double *dXdt);
void GetCodedJacobian(double t, const double *X, double *dfdX, double *dfdt);
void PlotVelocityFields(double t);
void Setup_PlotVelocity(long int *pktdetails);
void CleanUp_PlotVelocity();
void CombinePlotFiles(int numpackets);
void GeneratePlotMesh(long int *pktdetails);
void CopyPlottoWork(long int pktnum);




//find index by bisection
inline int FindSlot(double x,int d, int *dres)  {
		
  unsigned long jl,jm,ju;
  jl=0;
  ju=dres[d]-1;
  while (ju-jl > (unsigned long)1) {
    jm = (ju+jl) >> 1;
    if (x >= Data_Grid[d][jm]) jl=jm;
    else ju=jm;
  }
  
  return int(jl);
}


inline int FindSlot(double x,int d, double *loc, int *dres)  {
		
  unsigned long jl,jm,ju;
  jl=0;
  ju=dres[d]-1;
  while (ju-jl > (unsigned long)1) {
    jm = (ju+jl) >> 1;
    if (x >= Data_Grid[d][jm]) jl=jm;
    else ju=jm;
  }
  
  loc[d] = (x-Data_Grid[d][jl])/(Data_Grid[d][jl+1]-Data_Grid[d][jl]);
  
  return int(jl);
}


inline void GetIJloc(double *X, double *loc, int *ij) { 
  
  double xmod;
  for(int d=0;d<ND;++d) { 
    if(Data_Periodic[d]) {
      xmod = fmod( X[d]-Data_Min[d], Data_Period[d]);
      if(xmod<0) {
        xmod+=Data_Period[d];
      }
      xmod += Data_Min[d];
    }
    else {
      xmod = X[d];
    }
    
    if(Data_NonUniformGrid[d]) {
      ij[d] = FindSlot(xmod,d,loc,Data_Res);
    }  
    else {
      ij[d] = (int)floor((xmod-Data_Min[d])/Data_Delta[d]);
      if(ij[d] < 0) ij[d]=0;
      
      if(ij[d]>=(Data_Res[d]-1)){
        ij[d] =Data_Res[d]-2;
        loc[d] = 1.0;
      }
      else {
        loc[d] = ((xmod-Data_Min[d])-(ij[d]*Data_Delta[d]))/Data_Delta[d];
      }
    }
  }
}

inline void GetNestedIJloc(double *X, double *loc, int *ij, int nest, int tf) { 

  for(int d=0;d<ND;++d) { 
    
    if(Data_NonUniformGrid[d]) {
      ij[d] = FindSlot(X[d],d,loc,&DataNest_Res[tf][nest*ND]);
    }  
    else {
      ij[d] = (int)floor((X[d]-DataNest_Min[tf][nest*ND+d])/DataNest_Delta[tf][nest*ND+d]);
      if(ij[d] < 0) ij[d]=0;
      
      if(ij[d]>=(DataNest_Res[tf][nest*ND+d]-1)){
        ij[d] = DataNest_Res[tf][nest*ND+d]-2;
        loc[d] = 1.0;
      }
      else {
        loc[d] = ((X[d]-DataNest_Min[tf][nest*ND+d])-(ij[d]*DataNest_Delta[tf][nest*ND+d]))/DataNest_Delta[tf][nest*ND+d];
      }
    }
  }
}




inline void GetVelocity(double t,double *X,double *dXdt)  {
  
  switch (Velocity_Format) {
    case 0:
      GetCodedVelocity(t,X,dXdt);
      return;
    case 1:
      GetEquationVelocity(t,X,dXdt);
      return;
    case 2:
      GetCartesianVelocity(t,X,dXdt);
      return;
    default:
      GetCodedVelocity(t,X,dXdt);
      return;
  }
}


#endif


