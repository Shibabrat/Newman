// ftle.h
//
// Header file for the ftle module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _FTLE_H_
#define _FTLE_H_

#include <cmath>
#include "macros.h"
#include "indices.h"
#include <iostream>
#include <map>
#include <vector>
using namespace std;

extern double FBP_mu;

extern  int   myrank;
extern int    MapCoord;

extern int    LCS_Extract;
extern int    LCS_NumFields;

extern int    Track_Storm;
extern char   Track_Input[LONGSTRING]; 

extern int    Velocity_Format;
extern double Data_TMin;
extern double Data_TMax;
extern int    Data_TPeriodic;
extern int    Data_TRes;

extern int    Output_ForeignEndian;
extern double Output_T1;
extern double Output_TDelta;
extern int    Output_TRes;
extern int    MyOutput_TRes;

extern int    Atmos_Set;
extern double Atmos_Radius;
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
extern int FTLE_Type;

extern char FTLE_OutFile[SHORTSTRING];
extern char FTLE_Output[LONGSTRING];
extern char FTLE_Work[LONGSTRING];
extern char FTLE_Scratch[LONGSTRING];

extern char Path_Output[LONGSTRING];

struct FTLE_Point {
  double FTLEValue;
  int    HaveFTLE;
};

struct FTLE_Drift {
  double X[ND];
};

struct FTLE_Slide {
  int    slidenumber;
  int    blocksize;
  int    numdrifters;
  double launchtime;
  double stoptime;
  double completestatus;
  double ftlemin[ND];
};

extern struct FTLE_Slide *FTLE;

extern struct FTLE_Point *FTLE_Pts;

extern map< int, struct FTLE_Drift> FTLE_Dfts;

extern int    FTLE_BlockSize;


// Function Prototypes:
void CombineFTLEFiles(int ss);
void OutputFTLEtoFile(int ss);
int  ReadInTrack(vector<double>& trackarray);
void GetFTLE(int index, double FTLE_IntTime, int j1, int j2);
void AllocateMemoryForFTLE(long int *pktdetails);
void FreeMemoryForFTLE();
double GetFTLEForPoint(int index, double FTLEtime);
void ComputeFTLE(int ss, long int j1, long int packetID);
void ReadInFTLESlide(int ss, long int packetID);
void WriteOutFTLESlide(int ss, long int packetID);
int Get_OmegaEval(double FTLEarray[], double omega[], double mineval[], double smFTLE[]);
void CopyFTLEtoWork(long int pktnum);

inline void BuildDifferenceMatrix(int f0, double pts[ND][2][ND]) {
  
  int ij0[ND];
  int ijn[ND];
  
  f2ij(f0,ij0,FTLE_Res); 
  
  for(int d=0;d<ND;++d) {  
    
    for(int di=0;di<ND;++di)
      ijn[di] = ij0[di]+1;
    
    ijn[d] = ij0[d];
     
    for(int di=0;di<ND;++di) 
      pts[d][0][di] = FTLE_Dfts[ij2f(ijn,FTLE_DftRes)].X[di];
    
    ijn[d] = ij0[d] + 2;
    
    for(int di=0;di<ND;++di)
      pts[d][1][di] = FTLE_Dfts[ij2f(ijn,FTLE_DftRes)].X[di];
    
   }      
}

  

#endif


