// globals.h
//
//  Declaration of global variables
//  Global variables are deprecated, but what I am I supposed to do?
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "macros.h"
#include <fstream>

vector<int> Time_Origin;

double FBP_mu;
double FBP_a;
double FBP_omega;
double FBP_theta0;
double FBP_mu3;
double FBP_sunbdy2;
double FBP_moonbdy2;
double FBP_earthbdy2;

double Eye_Amplitude;
double Eye_Omega;
double Eye_A;
double Eye_alpha;
double Eye_y0;

double  ABC_Amplitude;

double Pendulum_Amplitude;

double **Utransform;

int     DNA_N;
vector<double>  DNA_q0;

double Velocity_Null[ND];

int    myrank;
int    NumProcs;

int    Query_Velocity;
double Query_X[ND];

double Parallel_LoadRatio;

int    Nest_NumNests;
vector<int> Nest_List;

int    MapCoord;
double TBP_mu;
double TBP_ecc;
int     TBP_Ordinate;
int     TBP_VDirection;

int LCS_Extract;
int LCS_NumFields;

struct FTLE_Slide *FTLE;
struct FTLE_Point *FTLE_Pts;
map<int,struct FTLE_Drift> FTLE_Dfts;

double **Data_Array;
double ***DataNest_Array;
double **DataNest_Min;
double **DataNest_Max;
double **DataNest_Delta;
int   **DataNest_Res;

int    Track_Storm;
double *Track_Array;

double Boundary_MaskValue;
int    Boundary_Method;
int    Boundary_NumLists;

int    Velocity_Format;
int    Data_NumInputFiles;
int    Data_Format;
double Data_Min[ND];
double Data_Max[ND];
int    Data_Res[ND];
double Data_Delta[ND];
int    Data_Periodic[ND];
double Data_Period[ND];
int    Data_NonUniformGrid[ND];
double *Data_Grid[ND];
int    Data_BlockSize;
int    Data_HeaderSize;

double Data_TMin;
double Data_TMax;
int    Data_TRes;
double Data_TDelta;
int    Data_TPeriodic;
double Data_TPeriod;

double datatime1;
double datatime2;

int    Data_FirstFrame;
int    Data_LastFrame;

int    Cell_NumVerts;
double *Cell_Areas;
int    Cell_Mask[ND];

int    Atmos_Set;
double Atmos_Radius;

int    FTLE_Compute;
int    Trace_Compute;

double Int_T1;
double Int_T2;
double Output_T1;
double Output_T2; 
double Output_TDelta;
int    Output_TRes;
int    MyOutput_TRes;
double nextoutputtime;

double FTLE_IntTLength;
int    Time_Direction;


int    FTLE_Type;
double FTLE_TrackWidth[ND];
double FTLE_Min[ND];
double FTLE_Max[ND];
int    FTLE_Res[ND];
int    FTLE_DftRes[ND];
double FTLE_Delta[ND];
double FTLE_Offset;
int    FTLE_BlockSize;
double Filter_Width;

int    Int_Method;
double Int_AbsTol;
double Int_RelTol;
double Int_TimeStep;
double Int_MinTimeStep;
double Int_MaxTimeStep;
int    Int_UseJacobian;

int    Trace_GenerateMesh;
int    Trace_NumDrifters;
double Trace_MeshMin[ND]; 
double Trace_MeshMax[ND];
int    Trace_MeshRes[ND];
double Trace_MeshDelta[ND];

double Trace_MeshReleaseTime;
double Trace_MeshColor;
int    Trace_ColorDimension;

int    Plot_Velocity;
double Plot_Min[ND];
double Plot_Max[ND];
int    Plot_Res[ND];
double Plot_Delta[ND];
int    Plot_BlockSize;

int    Parameters_Print;

char AnalyticEq[ND][LONGSTRING];

char FTLE_OutFile[SHORTSTRING];
char FTLE_Output[LONGSTRING];
char FTLE_Work[LONGSTRING];
char FTLE_Scratch[LONGSTRING];
char Data_InFile[SHORTSTRING];
char Data_Input[LONGSTRING];
char Data_Work[LONGSTRING];
char Data_Scratch[LONGSTRING];
char Trace_InFile[SHORTSTRING];
char Trace_Input[LONGSTRING];
char Trace_Work[LONGSTRING];
char Trace_OutFile[SHORTSTRING];
char Trace_Output[LONGSTRING];
char Trace_Scratch[LONGSTRING];
char Boundary_InFile[SHORTSTRING];
char Boundary_Input[SHORTSTRING];
char Track_InFile[SHORTSTRING];
char Track_Input[LONGSTRING]; 
char Plot_OutFile[SHORTSTRING];
char Plot_Output[LONGSTRING];
char Plot_Work[LONGSTRING];
char Plot_Scratch[LONGSTRING];

char Path_Work[LONGSTRING];
char Path_Input[LONGSTRING];
char Path_Output[LONGSTRING];
char Path_Scratch[LONGSTRING];

#endif
