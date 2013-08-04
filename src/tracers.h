//tracers.h
//
// Header file for the tracers module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _TRACERS_H_
#define _TRACERS_H_


#include <iostream>
#include <fstream>
using namespace std;
#include "boundary.h"
#include "integrate.h"
#include "indices.h"

extern vector<int> Time_Origin;

extern int    myrank;
extern int    Output_TRes;

extern int    Atmos_Set;
extern double Atmos_Radius;

extern int    Trace_GenerateMesh;
extern int    Trace_NumDrifters;

extern double Trace_MeshMin[ND]; 
extern double Trace_MeshMax[ND];
extern int    Trace_MeshRes[ND];
extern double Trace_MeshDelta[ND];

extern double Trace_MeshReleaseTime;
extern double Trace_MeshColor;
extern int    Trace_ColorDimension;

extern char Trace_InFile[SHORTSTRING];
extern char Trace_Input[LONGSTRING];
extern char Trace_OutFile[SHORTSTRING];
extern char Trace_Output[LONGSTRING];
extern char Trace_Scratch[LONGSTRING];
extern char Trace_Work[LONGSTRING];

extern double Int_T1;
extern double Output_T1;
extern int  Time_Direction;

struct Trace_Point {
  double X[ND];
  double ReleaseTime;
  double Color;
  int    LeftDomain;
};

// Function Prototypes:
void CombineTraceFiles(int numprocs);

void OutputTracers(double t);
void OutputTracersTecplotASCII(double t);
void OutputTracersRawBinary(double t);
void IntegrateTracers(double t1, double t);
void ReadInTracerData(long int *pktdetails);
void GenerateTracerMesh(long int *pktdetails);
void Setup_Trace(long int *pktdetails);
void CleanUpTrace();
void CopyTracetoWork(long int pktnum);


#endif

