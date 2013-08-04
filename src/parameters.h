// parameters.h
//
// Header file for the parameters module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <cmath>
#include "macros.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
using namespace std;

extern vector<int> Time_Origin;

extern double FBP_mu;
extern double FBP_a;
extern double FBP_omega;
extern double FBP_theta0;
extern double FBP_mu3;
extern double FBP_sunbdy2;
extern double FBP_moonbdy2;
extern double FBP_earthbdy2;

extern double Eye_Amplitude;
extern double Eye_Omega;
extern double Eye_A;
extern double Eye_alpha;
extern double Eye_y0;

extern double Pendulum_Amplitude;
extern double  ABC_Amplitude;
extern double Velocity_Null[ND];

extern  int    Query_Velocity;
extern  double Query_X[ND];

extern vector<int> Nest_List;
extern  int   myrank;

extern int    Nest_NumNests;

extern double Parallel_LoadRatio;
extern int    LCS_Extract;
extern int    LCS_NumFields;

extern int    Cell_NumVerts;
extern double *Cell_Areas;
extern int    Cell_Mask[ND];

extern int    Parameters_Print;

extern int    Track_Storm;
extern double *Track_Array;

extern int    DNA_N;

extern int    MapCoord;
extern double TBP_mu;
extern double TBP_ecc;
extern int    TBP_Ordinate;
extern int    TBP_VDirection;

extern double Boundary_MaskValue;
extern int    Boundary_Method;
extern int    Boundary_NumLists;

extern int    Velocity_Format;
extern char AnalyticEq[ND][LONGSTRING];

extern int    Data_NumInputFiles;
extern int    Data_Format;

extern double Data_TMin;
extern double Data_TMax;
extern int    Data_TRes;
extern double Data_TDelta;

extern double Data_TMinReq;
extern double Data_TMaxReq;
extern int    Data_TPeriodic;
extern double Data_TPeriod;
extern int    Data_Periodic[ND];
extern double Data_Period[ND];
extern int    Data_NonUniformGrid[ND];

extern double datastart;
extern double dataend;
extern int    df;
extern int    Data_FirstFrame;
extern int    Data_LastFrame;

extern int    Atmos_Set;
extern double Atmos_Radius;

extern int    FTLE_Compute;
extern int    Trace_Compute;

extern double Output_T1;
extern double Output_T2; 
extern double Int_T1;
extern double Int_T2;
extern double Output_TDelta;
extern int    Output_TRes;
extern int    MyOutput_TRes;
extern double nextoutputtime;

extern double FTLE_IntTLength;
extern int    Time_Direction;

extern double FTLE_TrackWidth[ND];
extern double FTLE_Min[ND];
extern double FTLE_Max[ND];
extern int    FTLE_Res[ND];
extern int    FTLE_DftRes[ND];
extern double FTLE_Delta[ND];
extern double FTLE_Offset;
extern int    FTLE_BlockSize;
extern double Filter_Width;
extern int    FTLE_Type;

extern int    Int_Method;
extern double Int_AbsTol;
extern double Int_RelTol;
extern int    Int_UseJacobian;
extern double Int_TimeStep;
extern double Int_MinTimeStep;
extern double Int_MaxTimeStep;

extern int    Trace_GenerateMesh;
extern int    Trace_NumDrifters;
extern double Trace_MeshMin[ND]; 
extern double Trace_MeshMax[ND];
extern int    Trace_MeshRes[ND];
extern double Trace_MeshDelta[ND];

extern double Trace_MeshReleaseTime;
extern double Trace_MeshColor;
extern int    Trace_ColorDimension;

extern int    Plot_Velocity;
extern double Plot_Min[ND];
extern double Plot_Max[ND];
extern int    Plot_Res[ND];
extern double    Plot_Delta[ND];
extern int    Plot_BlockSize;

extern char Path_Work[LONGSTRING];
extern char Path_Input[LONGSTRING];
extern char Path_Output[LONGSTRING];
extern char Path_Scratch[LONGSTRING];

extern char FTLE_OutFile[SHORTSTRING];
extern char FTLE_Output[LONGSTRING];
extern char FTLE_Work[LONGSTRING];
extern char FTLE_Scratch[LONGSTRING];
extern char Data_InFile[SHORTSTRING];
extern char Data_Input[LONGSTRING];
extern char Data_Work[LONGSTRING];
extern char Data_Scratch[LONGSTRING];
extern char Trace_InFile[SHORTSTRING];
extern char Trace_Input[LONGSTRING];
extern char Trace_OutFile[SHORTSTRING];
extern char Trace_Output[LONGSTRING];
extern char Trace_Scratch[LONGSTRING];
extern char Boundary_InFile[SHORTSTRING];
extern char Boundary_Input[SHORTSTRING];
extern char Track_InFile[SHORTSTRING];
extern char Track_Input[LONGSTRING]; 
extern char Trace_Work[LONGSTRING];
extern char Plot_OutFile[SHORTSTRING];
extern char Plot_Output[LONGSTRING];
extern char Plot_Work[LONGSTRING];
extern char Plot_Scratch[LONGSTRING];

// Function Prototypes:
void SetDefaultParameters();
void ReadInParameters(char *argv[]);
void SetDerivedParameters();
void SetMyDerivedParameters(long int *packetdetails);
void CopyDatatoScratch(int datafirstframe,int datalastframe);
void CheckParameters();

#endif


