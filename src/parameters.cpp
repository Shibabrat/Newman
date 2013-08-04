// paraemters.cpp
//
// parameters module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//


#include "parameters.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void SetDefaultParameters() {
  
  if(myrank == 0) cout << "Setting default parameter values ... " << flush;
  
  for(int d=0;d<ND;++d)
    Velocity_Null[d] = 0;
  
  Pendulum_Amplitude = 0.0;
  DNA_N  = 30;
  Eye_y0 = -1.5;
  Eye_A = 0.716077;
  Eye_alpha = 0.15;
  Eye_Amplitude = 0.0;
  Eye_Omega = 1.0;
  ABC_Amplitude = 0.0;
  Query_Velocity = 0;
  Parameters_Print = 1;
  Parallel_LoadRatio = 1.0;
  FTLE_Compute = 1;
  FTLE_Type = 0;
  Trace_Compute = 0;
  Plot_Velocity = 0;
  Time_Direction = 1;
  Output_T1 = 0.0;
  Output_TRes = 1;
  MyOutput_TRes = 1;
  Output_TDelta = 1.0;
  Int_Method = 1;
  Int_AbsTol = 1e-9;
  Int_RelTol = 1e-6;
  Int_UseJacobian = 0;
  Int_TimeStep = 1e-4;
  Int_MinTimeStep = 0.0;
  Int_MaxTimeStep = 10.0;
  strcpy(Path_Input,"./");
  strcpy(Path_Output,"./");
  strcpy(Path_Work,"./");
  strcpy(Path_Scratch,"./");
  Velocity_Format = 0;
  for(int d=0;d<ND;++d)
    strcpy(AnalyticEq[d],"0.0");
  strcpy(Data_InFile,"velocitydata");
  Data_NumInputFiles = 0;
  Data_Format = 2;
  for(int d=0;d<ND;++d) {
    Data_Periodic[d] = 0;
    Data_NonUniformGrid[d] = 0;
  }
  Data_TRes = 2;
  Data_TDelta = 1.0;
  Data_TMin = 0.0;
  Data_TPeriodic = 0;  
  Atmos_Set = 0;
  Atmos_Radius = 6378.1;
  for(int d=0;d<ND;++d) {
    FTLE_Min[d] = 0.0;
    FTLE_Max[d] = 1.0;
    FTLE_Res[d] = 1;
  }
  MapCoord = 0;
  FTLE_IntTLength = 1.0;
  strcpy(FTLE_OutFile,"repatt");
  Filter_Width = 1.0;
  Trace_GenerateMesh = 0;
  strcpy(Trace_InFile,"drifterinput.dat");
  strcpy(Trace_OutFile,"drifterout");
  for(int d=0;d<ND;++d) {
    Trace_MeshMin[d] = 0.0; 
    Trace_MeshMax[d] = 1.0;
    Trace_MeshRes[d] = 2;
  }
  Trace_MeshReleaseTime = 0.0;
  Trace_MeshColor = 1;
  Trace_ColorDimension = (ND-1);
  strcpy(Plot_OutFile,"velocityplot");
  for(int d=0;d<ND;++d) {
    Plot_Min[d]=0.0;
    Plot_Max[d]=1.0;
    Plot_Res[d]=2;
  }
  Boundary_Method = 0;
  Boundary_MaskValue = -9999.0;
  strcpy(Boundary_InFile,"boundary.dat");
  Track_Storm = 0;
  strcpy(Track_InFile,"track.dat");
  TBP_mu = 0.1;
  TBP_ecc = 0.04;
  TBP_Ordinate = 0;
  TBP_VDirection = 1;
  LCS_Extract = 0;
  FBP_mu = 3.036e-6;
  FBP_a = 2.57356507e-3;
  FBP_omega = 12.3688695;
  FBP_theta0 = 0.0;
  FBP_mu3 = 3.733999e-8;
  FBP_sunbdy2 = 0.0;
  FBP_moonbdy2 = 0.0;
  FBP_earthbdy2 = 0.0;
 
  
  for(int d=0;d<6;++d)
    Time_Origin.push_back(0);
  
  if(myrank == 0) cout << "OK!" << endl;
  
}

void ReadInParameters(char *argv[]) {
  
  char buf[LONGSTRING];
  char pname[LONGSTRING];
  
  ifstream fp(argv[1]);
  if(!fp.is_open()) {
    cout << "\n Error opening file: \a" << argv[1] << endl;
    exit(FAILURE);
  }
  
  if(myrank == 0) {
    cout << "Reading in Parameters from " << argv[1] << " ..." << endl;
  }
  
#define pget(a,b) if(!strcmp(a,buf)) { fp >> (b); continue; }
#define pgetld(a,b) if(!strcmp(a,buf)) { fp.getline(b,LONGSTRING); dflag=1; break; }
#define pgetl(a,b) if(!strcmp(a,buf)) { fp.getline(b,LONGSTRING); continue; }
#define pgetd(a,b) if(!strcmp(a,buf)) { fp >> (b);  dflag=1; break;  }
  
  int dflag=0;
  while( (fp >> buf) != NULL ) {
    
    if(buf[0] == '#') {
      fp.ignore(LONGSTRING,'\n');
      continue;
    }
    pget("FBP_mu",FBP_mu);
    pget("FBP_mu3",FBP_mu3);
    pget("FBP_a",FBP_a);
    pget("FBP_omega",FBP_omega);
    pget("FBP_theta0",FBP_theta0);
    pget("FBP_sunbdy2",FBP_sunbdy2);
    pget("FBP_earthbdy2",FBP_earthbdy2);
    pget("FBP_moonbdy2",FBP_moonbdy2);
    pget("Pendulum_Amplitude",Pendulum_Amplitude);
    pget("DNA_N",DNA_N);
    pget("Eye_alpha",Eye_alpha);
    pget("Eye_A",Eye_A);
    pget("Eye_y0",Eye_y0);
    pget("Eye_Omega",Eye_Omega);
    pget("Eye_Amplitude",Eye_Amplitude);
    pget("Atmos_Radius",Atmos_Radius);
    pget("Atmos_Set",Atmos_Set);
    pget("Boundary_InFile",Boundary_InFile);
    pget("Boundary_Method",Boundary_Method);
    pget("Boundary_MaskValue",Boundary_MaskValue);
    pget("Data_Format",Data_Format);
    pget("Data_InFile",Data_InFile);
    for(int d=0;d<ND;++d) {
      sprintf(pname,"Data_Periodic[%d]",d);
      pgetd(pname,Data_Periodic[d]);
      sprintf(pname,"Data_NonUniformGrid[%d]",d);
      pgetd(pname,Data_NonUniformGrid[d]);
    }
    if(dflag) {dflag=0;continue;}
    
     pget("Data_NumInputFiles",Data_NumInputFiles);
    if(!strcmp("Time_Origin",buf)) {
      Time_Origin.clear();
      char *nptr;
      char *eptr;
      fp.getline(buf,LONGSTRING);
      nptr = &buf[0];
      int tval = (int)strtod(nptr,&eptr);
      while(*nptr != '\0') {
        nptr=eptr;
        Time_Origin.push_back(tval);
        tval = (int)strtod(nptr,&eptr);
      }
      if(Time_Origin.size()!=6) {
        cout << "Time_Origin requires 6 int arguments: Year Month Day Hour Minute Second" << endl;
        exit(1);
      }
      continue;
    }
    
    pget("Data_TRes",Data_TRes);
    pget("Data_TDelta",Data_TDelta);
    pget("Data_TMin",Data_TMin);
    pget("Data_TPeriodic",Data_TPeriodic);
    pget("Filter_Width",Filter_Width); 
    pget("FTLE_Compute",FTLE_Compute);
    pget("FTLE_Type",FTLE_Type);
    for(int d=0;d<ND;++d) {
      sprintf(pname,"FTLE_Min[%d]",d);
      pgetd(pname,FTLE_Min[d]);
      sprintf(pname,"FTLE_Max[%d]",d);
      pgetd(pname,FTLE_Max[d]);
      sprintf(pname,"FTLE_Res[%d]",d);
      pgetd(pname,FTLE_Res[d]);
    }
    if(dflag) {dflag=0;continue;}
     for(int d=0;d<ND;++d) {
      sprintf(pname,"FTLE_TrackWidth[%d]",d);
      pgetd(pname,FTLE_TrackWidth[d]);
    }
    if(dflag) {dflag=0;continue;}
    pget("FTLE_IntTLength",FTLE_IntTLength);
    pget("FTLE_OutFile",FTLE_OutFile);

    pget("Int_Method",Int_Method);
    pget("Int_AbsTol",Int_AbsTol);
    pget("Int_RelTol",Int_RelTol);
    pget("Int_UseJacobian",Int_UseJacobian);
    pget("Int_TimeStep",Int_TimeStep);
    pget("Int_MinTimeStep",Int_MinTimeStep);
    pget("Int_MaxTimeStep",Int_MaxTimeStep);
    
    pget("LCS_Extract",LCS_Extract);
    pget("MapCoord",MapCoord);
    
    if(!strcmp("Nest_List",buf)) {
      char *nptr;
      char *eptr;
      fp.getline(buf,LONGSTRING);
      nptr = &buf[0];
      int nestnum = (int)strtod(nptr,&eptr);
      while(*nptr != '\0') {
        nptr=eptr;
        Nest_List.push_back(nestnum);
        nestnum = (int)strtod(nptr,&eptr);
      }
      continue;
    }
    pget("Nest_NumNests",Nest_NumNests);

    pget("Output_T1",Output_T1);
    pget("Output_TRes",Output_TRes);
    pget("Output_TDelta",Output_TDelta);
    
    pget("Parallel_LoadRatio",Parallel_LoadRatio);
    
    pget("Parameters_Print",Parameters_Print);
    pget("Path_Input",Path_Input);
    pget("Path_Output",Path_Output);
    pget("Path_Work",Path_Work);
    pget("Path_Scratch",Path_Scratch);
    for(int d=0;d<ND;++d) {
      sprintf(pname,"Plot_Min[%d]",d);
      pgetd(pname,Plot_Min[d]);
      sprintf(pname,"Plot_Max[%d]",d);
      pgetd(pname,Plot_Max[d]);
      sprintf(pname,"Plot_Res[%d]",d);
      pgetd(pname,Plot_Res[d]);
    }
    if(dflag) {dflag=0;continue;}
    pget("Plot_Velocity",Plot_Velocity);
    pget("Plot_OutFile",Plot_OutFile);
    pget("Query_Velocity",Query_Velocity);
    for(int d=0;d<ND;++d) {
      sprintf(pname,"Query_X[%d]",d);
      pgetd(pname,Query_X[d]);
    }
    if(dflag) {dflag=0;continue;} 
    
    pget("Time_Direction",Time_Direction);
    pget("Trace_Compute",Trace_Compute);
    pget("Trace_GenerateMesh",Trace_GenerateMesh);
    pget("Trace_InFile",Trace_InFile);
    pget("Trace_OutFile",Trace_OutFile);
    for(int d=0;d<ND;++d) {
      sprintf(pname,"Trace_MeshMin[%d]",d);
      pgetd(pname,Trace_MeshMin[d]);
      sprintf(pname,"Trace_MeshMax[%d]",d);
      pgetd(pname,Trace_MeshMax[d]);
      sprintf(pname,"Trace_MeshRes[%d]",d);
      pgetd(pname,Trace_MeshRes[d]);
    }
    if(dflag) {dflag=0;continue;}
    pget("Trace_MeshReleaseTime",Trace_MeshReleaseTime);
    pget("Trace_MeshColor",Trace_MeshColor);
    pget("Trace_ColorDimension",Trace_ColorDimension);
    
    pget("Track_Storm",Track_Storm); 
    pget("Track_InFile",Track_InFile);
    
    for(int d=0;d<ND;++d) {
      sprintf(pname,"V[%d]",d);
      pgetld(pname,AnalyticEq[d]);
    }
    if(dflag) {dflag=0;continue;} 
     pget("Velocity_Format",Velocity_Format);
    
    for(int d=0;d<ND;++d) {
      sprintf(pname,"Velocity_Null[%d]",d);
      pgetd(pname,Velocity_Null[d]);
    }
    if(dflag) {dflag=0;continue;}
    
    pget("ABC_Amplitude",ABC_Amplitude);
    pget("TBP_mu",TBP_mu);
    pget("TBP_ecc",TBP_ecc);
    pget("TBP_Ordinate",TBP_Ordinate);
    pget("TBP_VDirection",TBP_VDirection);
   
    if(myrank==0) cout << "Unknown parameter: " << buf << endl;
    exit(1);
  }
  
  fp.close();
  
#define pout(a,b) cout << a << "\t\t\t " << b << endl;
#define poutd(a,b) cout << a << "[" << d << "]" << "\t\t\t " << b << endl;
  
  if((myrank == 0) && Parameters_Print){
    pout("FBP_mu",FBP_mu);
    pout("FBP_mu3",FBP_mu3);
    pout("FBP_a",FBP_a);
    pout("FBP_omega",FBP_omega);
    pout("FBP_theta0",FBP_theta0);
    pout("FBP_sunbdy2",FBP_sunbdy2);
    pout("FBP_earthbdy2",FBP_earthbdy2);
    pout("FBP_moonbdy2",FBP_moonbdy2);
    pout("Pendulum_Amplitude",Pendulum_Amplitude);
    pout("DNA_N",DNA_N);
    pout("Eye_alpha",Eye_alpha);
    pout("Eye_A",Eye_A);
    pout("Eye_y0",Eye_y0);
    pout("Eye_Omega",Eye_Omega);
    pout("Eye_Amplitude",Eye_Amplitude);
    pout("Nest_NumNests",Nest_NumNests);
    pout("Parallel_LoadRatio",Parallel_LoadRatio);
    pout("LCS_Extract",LCS_Extract);
    pout("Path_Input",Path_Input);
    pout("Path_Output",Path_Output);
    pout("Path_Work",Path_Work);
    pout("Path_Scratch",Path_Scratch);
    pout("FTLE_Compute",FTLE_Compute);
    pout("FTLE_Type",FTLE_Type);
    pout("Trace_Compute",Trace_Compute);
    pout("Plot_Velocity",Plot_Velocity);
    pout("Velocity_Format",Velocity_Format);
    for(int d=0;d<ND;++d)
      poutd("AnalyticEq",AnalyticEq[d]);
    pout("Data_InFile",Data_InFile);
    pout("Data_NumInputFiles",Data_NumInputFiles);
    pout("Data_Format",Data_Format);
    for(int d=0;d<ND;++d) {
      poutd("Data_Periodic",Data_Periodic[d]);
      poutd("Data_NonUniformGrid",Data_NonUniformGrid[d]);
    }
    for(int d=0;d<ND;++d)
      poutd("Velocity_Null",Velocity_Null[d]);
    pout("Data_TRes",Data_TRes);
    pout("Data_TDelta",Data_TDelta);
    pout("Data_TMin",Data_TMin);
    pout("Data_TPeriodic",Data_TPeriodic);
    pout("Atmos_Set",Atmos_Set);
    pout("Atmos_Radius",Atmos_Radius);
    pout("Output_T1",Output_T1);
    pout("Output_TRes",Output_TRes);
    pout("Output_TDelta",Output_TDelta);
    pout("Int_Method",Int_Method);
    pout("Int_AbsTol",Int_AbsTol);
    pout("Int_RelTol",Int_RelTol);
    pout("Int_UseJacobian",Int_UseJacobian);
    pout("Int_TimeStep",Int_TimeStep);
    pout("Int_MinTimeStep",Int_MinTimeStep);
    pout("Int_MaxTimeStep",Int_MaxTimeStep);
    pout("Time_Direction",Time_Direction);
    for(int d=0;d<ND;++d) {
      poutd("FTLE_Min",FTLE_Min[d]);
      poutd("FTLE_Max",FTLE_Max[d]);
      poutd("FTLE_Res",FTLE_Res[d]);
      poutd("FTLE_TrackWidth",FTLE_TrackWidth[d]);
    }
    pout("FTLE_IntTLength",FTLE_IntTLength);
    pout("FTLE_OutFile",FTLE_OutFile);
    pout("Filter_Width",Filter_Width);
    pout("Trace_GenerateMesh",Trace_GenerateMesh);
    pout("Trace_InFile",Trace_InFile);
    pout("Trace_OutFile",Trace_OutFile);
    for(int d=0;d<ND;++d) {
      poutd("Trace_MeshMin",Trace_MeshMin[d]);
      poutd("Trace_MeshMax",Trace_MeshMax[d]);
      poutd("Trace_MeshRes",Trace_MeshRes[d]);
    }
    pout("Trace_MeshReleaseTime",Trace_MeshReleaseTime);
    pout("Trace_MeshColor",Trace_MeshColor);
    pout("Trace_ColorDimension",Trace_ColorDimension);
    pout("Plot_Velocity",Plot_Velocity);
    pout("Plot_OutFile",Plot_OutFile);
    for(int d=0;d<ND;++d) {
      poutd("Plot_Min",Plot_Min[d]);
      poutd("Plot_Max",Plot_Max[d]);
      poutd("Plot_Res",Plot_Res[d]);
    }
    pout("Boundary_Method",Boundary_Method);
    pout("Boundary_MaskValue",Boundary_MaskValue);
    pout("Track_Storm",Track_Storm);
    pout("Track_InFile",Track_InFile);
    pout("TBP_mu",TBP_mu);
    pout("TBP_ecc",TBP_ecc);
    pout("TBP_Ordinate",TBP_Ordinate);
    pout("TBP_VDirection",TBP_VDirection);
    pout("MapCoord",MapCoord);
    pout("Nest_NumNests",Nest_NumNests);
    cout << "Nest_List\t\t";
    vector<int>::iterator iter;
    for( iter = Nest_List.begin(); iter != Nest_List.end(); iter++)
      cout << *iter << " ";
    cout << endl;
    
    cout << "Time_Origin\t\t";
    for( iter = Time_Origin.begin(); iter != Time_Origin.end(); iter++)
      cout << *iter << " ";
    cout << endl;
  }
}

void SetMyDerivedParameters(long int *packetdetails) {
  
  if(myrank == 0) cout << "Setting local derived parameters ..." << flush;
  
  /*  Output_T1 is the time of the first output frame
  Output_T2 is the time of the last output frame
  Output_T2 >= Outout_T1 since we always view LCS movies playing forward in time
  */  
  int td = Time_Direction;  /* abbreviation for Time_Direction */
  
  int jbegin = packetdetails[1]; 
  int jend = packetdetails[2]; 
  
  if(FTLE_Compute) {
    int ssbegin = int(jbegin/FTLE_BlockSize);
    int ssend = int(jend/FTLE_BlockSize);
    
    MyOutput_TRes = ssend-ssbegin + 1;
    double MyOutput_T1, MyOutput_T2;
    MyOutput_T1 = Output_T1 + ssbegin*Output_TDelta;
    MyOutput_T2 = MyOutput_T1 + (MyOutput_TRes - 1) * Output_TDelta;
    Time_Direction > 0 ? (Int_T1 = MyOutput_T1) : (Int_T1 = MyOutput_T2);
    Time_Direction > 0 ? (Int_T2 = MyOutput_T2) : (Int_T2 = MyOutput_T1);
  }
  
  if(Trace_Compute || Plot_Velocity) {
    Int_T1 = Output_T1;
    Int_T2 = Output_T2;
  }
  
  nextoutputtime = Int_T1;
  
  if(FTLE_Compute) {
    
    if(Velocity_Format > 1) {
      double Data_T2Req;
      Data_T2Req = Int_T2 + Time_Direction * FTLE_IntTLength;
      if(!Data_TPeriodic) {
        if((Data_T2Req < Data_TMin)) Data_T2Req = Data_TMin;
        if((Data_T2Req > Data_TMax)) Data_T2Req = Data_TMax;
      }
      Data_LastFrame =  (int) (  td * floor( td * ( (Data_T2Req - Data_TMin) / Data_TDelta  - td * TINY_DATA ) ) );
      if(td*Data_LastFrame < td*Data_FirstFrame) Data_LastFrame = Data_FirstFrame;

    }
  }
  
  if(Velocity_Format > 1) 
    CopyDatatoScratch(Data_FirstFrame,Data_LastFrame);
  
  if(myrank == 0) cout << "OK!" << endl;
}


void CopyDatatoScratch(int datafirstframe,int datalastframe) {

  char syscommand[LONGSTRING];
  
  if(myrank == 0)
    cout << "\nCopying data to scratch file ... " << flush;
  int f0;
  int f1;
  if(datafirstframe < datalastframe) {
    f0 = datafirstframe; 
    f1 = datalastframe;
  }
  else {
    f0 = datalastframe; 
    f1 = datafirstframe;
  }
  for(int ff=f0;ff<=f1;++ff) {
    if(Velocity_Format==2) {
      sprintf(syscommand,"cp -u %sf%04d.new %s",Data_Work,ff,Path_Scratch);
      system(syscommand);
    }
    else if(Velocity_Format==3) { 
       vector<int>::iterator iter;
       for( iter = Nest_List.begin(); iter != Nest_List.end(); iter++) {
          int nn = *iter;
        sprintf(syscommand,"cp -u %sf%04d_%04d.new %s",Data_Work,ff,nn,Path_Scratch);
        system(syscommand);
      }
    }
  }
  if(myrank == 0)
    cout << "OK!\n" << flush;

}


void SetDerivedParameters() {
  
  /* Set path to various directories */
    if( Path_Input[strlen(Path_Input)-1] != '/' ) 
      sprintf(Path_Input,"%s/",Path_Input);
      
    if( Path_Work[strlen(Path_Work)-1] != '/' ) 
      sprintf(Path_Work,"%s/",Path_Work);
      
    if( Path_Output[strlen(Path_Output)-1] != '/' ) 
      sprintf(Path_Output,"%s/",Path_Output);
      
    if( Path_Scratch[strlen(Path_Scratch)-1] != '/' ) 
      sprintf(Path_Scratch,"%s/",Path_Scratch);
    
  /* Create various file strings */
    if(Velocity_Format > 1) {
      sprintf(Data_Input,"%s%s",Path_Input,Data_InFile);
      sprintf(Data_Work,"%s%s",Path_Work,Data_InFile);
      sprintf(Data_Scratch,"%s%s",Path_Work,Data_InFile);
    }
    
    if(Track_Storm)
      sprintf(Track_Input,"%s%s",Path_Input,Track_InFile);
    
    if(Trace_Compute) {
      sprintf(Trace_Input,"%s%s",Path_Input,Trace_InFile);
      sprintf(Trace_Scratch,"%strace",Path_Scratch);
      sprintf(Trace_Work,"%strace",Path_Work);
      sprintf(Trace_Output,"%s%s.raw",Path_Output,Trace_OutFile);
    }
    
    if(Plot_Velocity) {
      sprintf(Plot_Scratch,"%splot",Path_Scratch);
      sprintf(Plot_Work,"%splot",Path_Work);
      sprintf(Plot_Output,"%s%s.raw",Path_Output,Plot_OutFile);
    }
    
    if(Boundary_Method == 2)
      sprintf(Boundary_Input,"%s%s",Path_Input,Boundary_InFile);
 
    if(FTLE_Compute) {
      if(!strcmp(FTLE_OutFile,"repatt")) {
        (Time_Direction > 0) ? sprintf(FTLE_OutFile, "FTLErep") : sprintf(FTLE_OutFile, "FTLEatt");
      }
      sprintf(FTLE_Output,"%s%s",Path_Output,FTLE_OutFile);
      sprintf(FTLE_Work,"%snm",Path_Work);
      sprintf(FTLE_Scratch,"%snm",Path_Scratch);
    }
    
      
  /*  Output_T1 is the time of the first output frame
  Output_T2 is the time of the last output frame
  Output_T2 >= Outout_T1 since we always view LCS movies playing forward in time
  
  Int_T1 is the 
  */  
  int td = Time_Direction;  /* abbreviation for Time_Direction */
  
  if(myrank == 0) cout << "Setting derived parameters ..." << endl;
  
  if(FTLE_Compute) {
    Output_T2 = Output_T1 + (Output_TRes - 1) * Output_TDelta;
    Time_Direction > 0 ? (Int_T1 = Output_T1) : (Int_T1 = Output_T2);
    Time_Direction > 0 ? (Int_T2 = Output_T2) : (Int_T2 = Output_T1);
  }
  if(Trace_Compute || Plot_Velocity) {
    Output_T2 = Output_T1 + Time_Direction * (Output_TRes - 1) * Output_TDelta;
    Int_T1 = Output_T1;
    Int_T2 = Output_T2;
  }
  
  nextoutputtime = Int_T1;
  
  if(Velocity_Format > 1) {         
    if(Data_TDelta <= 0) {
      if(Data_TRes == 1) {
        Data_TDelta = 1;
      }
      else {
        cout << "Parameter Error: Data_TDelta <= 0\a" << endl;
        exit(1);
      }
    }
    
    Data_TMax = Data_TMin + (Data_TRes - 1) * Data_TDelta;
    Data_TPeriod = Data_TMax - Data_TMin;
    
    
    /* Set first data frame.  Next three lines must appear before the if(FTLE_Compute) conditional */
    Data_FirstFrame = (int) ( td * floor( td * ( (Int_T1 - Data_TMin) / Data_TDelta  + td * TINY_DATA ) ) );
    Data_LastFrame =  (int) ( td * floor( td * ( (Int_T2 - Data_TMin) / Data_TDelta  - td * TINY_DATA ) ) );
    if(td*Data_LastFrame < td*Data_FirstFrame) Data_LastFrame = Data_FirstFrame;
    
    Cell_NumVerts = 1;
    Cell_NumVerts <<= ND;
    Cell_Areas = new double[Cell_NumVerts];
    for(int d=0;d<ND;++d) {
      Cell_Mask[d] = 1;
      Cell_Mask[d] <<= d;
    }
  }
  
  if(FTLE_Compute) {
    
    LCS_Extract ? LCS_NumFields = 4 : LCS_NumFields = 1;
    
    FTLE_BlockSize = 1;
    for(int d=0;d<ND;++d)
      FTLE_BlockSize *= FTLE_Res[d];
    
    for(int d=0;d<ND;++d)
      FTLE_DftRes[d] = FTLE_Res[d]+2;
    
    for(int d=0;d<ND;++d) {
      if(Track_Storm && d<2) {
        if(FTLE_Res[d] == 1)
          FTLE_Delta[d] = FTLE_TrackWidth[d];
        else
          FTLE_Delta[d] = FTLE_TrackWidth[d]/(FTLE_Res[d]-1);
      }
      else {
        if(FTLE_Res[d] == 1)
          FTLE_Delta[d] = (FTLE_Max[d] - FTLE_Min[d]);
        else
          FTLE_Delta[d] = (FTLE_Max[d] - FTLE_Min[d]) / (FTLE_Res[d] - 1);
      }
    }
    
    if(Velocity_Format > 1) {
      double Data_T2Req;
      Data_T2Req = Int_T2 + Time_Direction * FTLE_IntTLength;
      if(!Data_TPeriodic) {
        if((Data_T2Req < Data_TMin)) Data_T2Req = Data_TMin;
        if((Data_T2Req > Data_TMax)) Data_T2Req = Data_TMax;
      }
      Data_LastFrame =  (int) ( td * floor( td * ( (Data_T2Req - Data_TMin) / Data_TDelta  - td * TINY_DATA ) ) );
      if(td*Data_LastFrame < td*Data_FirstFrame) Data_LastFrame = Data_FirstFrame;
    }
  }
  
  if(Trace_Compute) {
    
    if(Trace_GenerateMesh) {
      
      Trace_NumDrifters = 1;
      for(int d=0;d<ND;++d)
        Trace_NumDrifters *= Trace_MeshRes[d];
      
      for(int d=0;d<ND;++d) {
        if(Trace_MeshRes[d] == 1)
          Trace_MeshDelta[d] = 0.0;
        else
          Trace_MeshDelta[d] = (Trace_MeshMax[d] - Trace_MeshMin[d]) / ((double)Trace_MeshRes[d] - 1.0);
      }
    }
    else {
      
      ifstream Trace_InFileID(Trace_Input);
      if(!Trace_InFileID.is_open()) {
        cout << "Error opening file: \a" << Trace_Input << endl;
        exit(FAILURE);
      }
      
      char buf[LONGSTRING];
      Trace_InFileID.getline(buf,LONGSTRING,'\n');
      sscanf(buf, "%d\n", &Trace_NumDrifters);
      Trace_InFileID.close();
    }
  }
  
  if(Plot_Velocity) {
    Plot_BlockSize = 1;
    for(int d=0;d<ND;++d) {
      Plot_BlockSize *= Plot_Res[d];
      if(Plot_Res[d] == 1)
        Plot_Delta[d] = (Plot_Max[d]-Plot_Min[d]);
      else
        Plot_Delta[d] = (Plot_Max[d]-Plot_Min[d])/(Plot_Res[d] - 1);
    }
    
  }
  
  
  
  if(myrank == 0) cout << "OK!" << endl;
}

void CheckParameters() {
  
  if(myrank == 0) cout << "Checking parameters ... ";
  
  if(Output_TRes < 1) {
    cout << "Parameter Error: Output_TRes < 1. Nothing to do.\a" << endl;
    exit(1);
  }
  
  
  if(!Trace_Compute && !FTLE_Compute && !Plot_Velocity) {
    if(myrank == 0) cout << "Parameter Error: FTLE_Compute = Trace_Compute = Plot_Velocity = 0\a" << endl;
    exit(1);
  }
  
  if(Velocity_Format > 1) {
    
    if(!Data_TPeriodic) {
      if((Output_T1 < Data_TMin) || (Output_T1 > Data_TMax)) {
        cout << "Parameter Error: Output_T1 out of range\a" << endl;
        cout << Data_TMin << " : " << Output_T1 << " : " << Data_TMax << endl;
        exit(1);
      }
      if((Output_T2 < Data_TMin) || (Output_T2 > Data_TMax)) {
        cout << "Parameter Error: Output_T2 out of range\a" << endl;
        cout << Data_TMin << " : " << Output_T2 << " : " << Data_TMax << endl;
        exit(1);
      }
    }
    
  }
  
  if(Velocity_Format == 1)
    if(myrank == 0) cout << "\n\nATTENTION: Computation speed can be improved by coding the analytical equation directly in velocity.cpp\n" << endl;
  
  
  if(FTLE_Compute) {
    
    if( Output_T1 > Output_T2 ) {
      cout << "Parameter Error: Output_T1 > Output_T2\a" << endl;
      exit(1);
    }
    
    
    for(int d=0;d<ND;++d) {
      if(FTLE_Min[d] >= FTLE_Max[d]) {
        cout << "Parameter Error: FTLE_Min[" << d << "] >= FTLE_Max[" << d << "]\a" << endl;
        exit(1);
      }
      
      if(FTLE_Res[d] < 1) {
        cout << "Parameter Error: FTLE_Res[" << d << "] < 1\a" << endl;
        exit(1);
      }
    }
    if(!Data_TPeriodic && (Velocity_Format > 1)) {  
      if((Int_T2 + Time_Direction * FTLE_IntTLength) > Data_TMax || (Int_T2 + Time_Direction * FTLE_IntTLength) < Data_TMin) {
        cout << "Parameter WARNING: Insufficient Data to compute FTLE to time T for all slides!!\a" << endl;
      }
    }
  }
  
  if(Trace_Compute) {
    if(Trace_GenerateMesh) {
      if((Trace_MeshReleaseTime - Int_T1) * Time_Direction < 0) {
        cout << "\n\nParameter WARNING: Trace_MeshResleaseTime outside Output Times, setting to Int_T1.\a" << endl;
        Trace_MeshReleaseTime = Int_T1;
      }
      
      for(int d=0;d<ND;++d)
        if(Trace_MeshMin[d] > Trace_MeshMax[d])  
          cout << "Parameter WARNING: Trace_MeshMin[" << d << "] > Trace_MeshMax[" << d << "].\a" << endl;
    }
  }
  
  if(myrank == 0) cout << "OK!" << endl;
}
