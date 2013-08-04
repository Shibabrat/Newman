// data.h
//
// Header file for data module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#ifndef _DATA_H_
#define _DATA_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
using namespace std;
#include "macros.h"
#include "indices.h"

extern vector<int> Nest_List;

/* Globals that I need */
extern int    myrank;
extern int Nest_NumNests;
extern int Velocity_Format;

/* Parameters that I need from other modules */


/* Public parameters (other modules may need these parameters that I own) */


/* Private parameters (only I need these parameters) */
extern int    Data_NumInputFiles;
extern int    Data_Format;
extern int    Data_NumGrids;

extern double *Data_Grid[ND];
extern double Data_Min[ND];
extern double Data_Max[ND];
extern int    Data_Res[ND];
extern double Data_Delta[ND];
extern double Data_Period[ND];
extern int    Data_NonUniformGrid[ND];
extern int    Data_Periodic[ND];
extern int    Data_BlockSize;
extern int    Data_HeaderSize;

extern double Data_TMin;
extern double Data_TMax;
extern int    Data_TRes;
extern double Data_TDelta;

extern double **Data_Array;
extern double ***DataNest_Array;
extern double **DataNest_Min;
extern double **DataNest_Max;
extern double **DataNest_Delta;
extern int    **DataNest_Res;

extern double datastart;
extern double dataend;
extern int    Data_FirstFrame;
extern int    Data_LastFrame;

extern char Data_InFile[SHORTSTRING];
extern  char Data_Input[LONGSTRING];
extern  char Data_Work[LONGSTRING];
extern  char Data_Scratch[LONGSTRING];
extern  char Path_Input[LONGSTRING];

//Function prototypes
/*Public */

void LoadFirstDataFrame(int firstframe);
void UpdateDataFrame(int df2);
void UpdateDataNestFrame(int df2);
void LoadFirstDataNestFrame(int firstframe);

/* Private */
void ReadInNonUniformGrids();
void AllocateMemoryForVelocityData();
void FreeMemoryForVelocityData();
void ReadInASCIIVelocityData(int format, int numprocs);
void ReadInASCIIMultiLeanVelocityData(int nf);
void ReadInASCIIMultiLeanVelocityNestData(int nf, int Nest_NumNests);
void ReadInASCIITecplotVelocityData();
void ReadInDataParameters();
void ReadInDataNestParameters();

#endif


