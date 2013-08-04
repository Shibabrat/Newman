// tracers.cpp
//
// tracer module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include <cmath>
#include <vector>
#include "mapcoord.h"
#include "macros.h"
#include "tracers.h"
#include <fstream>

int    MyTrace_NumDrifters;
ofstream Trace_OutFileID;
struct Trace_Point *Trace_Array;

void CopyTracetoWork(long int pktnum) {

  char syscommand[LONGSTRING];
  sprintf(syscommand,"mv -f %s%04ld.raw %strace%04ld.raw",Trace_Scratch,pktnum,Path_Work,pktnum);
  system(syscommand);
    
}


void Setup_Trace(long int *pktdetails) {
  
  if(!Trace_GenerateMesh) 
    ReadInTracerData(pktdetails);
  else 
    GenerateTracerMesh(pktdetails);
  
  char tracefilename[LONGSTRING];
  //char syscommand[SHORTSTRING];
  //sprintf(syscommand,"rm -f %s*.raw",Trace_Scratch);
  //system(syscommand);
  sprintf(tracefilename,"%s%04ld.raw",Trace_Scratch,pktdetails[0]);
  Trace_OutFileID.open(tracefilename,ios::binary);
  if(!Trace_OutFileID.is_open()) {
    cout << "\n Error opening file: \a" <<  tracefilename << endl;
    exit(FAILURE);
  }
}

#include <sys/stat.h> 


void CombineTraceFiles(int numpackets) {
  
  char infilename[LONGSTRING];
  char systemcommand[LONGSTRING];
  
  vector<double> X;
  vector<double> color;
  int  numdrift;
  double frametime;
  
  if(myrank==0) cout << "Combining " << numpackets << " drifter files into raw binary format ... " << flush;
   
  ofstream fout(Trace_Output,ios::binary);
  if(!fout.is_open()) {
    cout << "\n Error opening file: \a" <<  Trace_Output << endl;
    exit(FAILURE);
  }

  int one = 1;
  fout.write((char *)&one, sizeof(int) );
  fout.write((char *)&Time_Origin[0], 6*sizeof(int) );
  fout.write((char *)&Output_TRes, sizeof(int) );
  fout.write((char *)&Atmos_Set, sizeof(int) );
  fout.write((char *)&Atmos_Radius, sizeof(double) );

  long int *fpos = new long int[numpackets];
  for(int p=0;p<numpackets;++p) fpos[p]=0;
  
  for(int tt=0;tt<Output_TRes;++tt) {  
    
    int count = 0;
    int colorcount = 0;
    int totnumdrift = 0;
    for(int d=0;d<ND;++d) {
      for(int p=0;p<numpackets;++p) {
        sprintf(infilename,"%s%04d.raw",Trace_Work,p);
        
        ifstream fin(infilename,ios::binary);
        if(!fin.is_open()) {
          cout << "\n Error opening file: \a" <<  infilename << endl;
          exit(FAILURE);
        }

        fin.seekg(fpos[p],ios::beg);

        fin.read((char *)&numdrift,sizeof(int));
        fin.read((char *)&frametime,sizeof(double));
        
        X.resize(count+numdrift,0);

        fin.read((char *)&X[count],numdrift * sizeof(double));
        
        count+=numdrift;
        if(d==0)
          totnumdrift+=numdrift;
        fpos[p]=fin.tellg();
        fin.close(); 
        
        
      }
    } 
    for(int p=0;p<numpackets;++p) {
      sprintf(infilename,"%s%04d.raw",Trace_Work,p);
      
      ifstream fin(infilename,ios::binary);
      if(!fin.is_open()) {
        cout << "\n Error opening file: \a" <<  infilename << endl;
        exit(FAILURE);
      }
      
      fin.seekg(fpos[p],ios::beg);
      
      fin.read((char *)&numdrift,sizeof(int));
      fin.read((char *)&frametime,sizeof(double));
      
      color.resize(colorcount+numdrift,0);
      
      fin.read((char *)&color[colorcount],numdrift * sizeof(double));
      colorcount+=numdrift;
      
      fpos[p]=fin.tellg();
      fin.close(); 
    }
    
    fout.write((char *)&colorcount,sizeof(int));
    fout.write((char *)&frametime,sizeof(double));
    fout.write((char *)&X[0],count*sizeof(double));
    fout.write((char *)&color[0],colorcount*sizeof(double));
    
  }
  
  fout.close();
  
  sprintf(systemcommand,"rm -f %s*.raw",Trace_Work);
  system(systemcommand);
  
  delete []fpos;
  cout << "\nOutput file: " << Trace_Output << endl; 
}


void CleanUpTrace() {
  
  delete []Trace_Array;
  Trace_OutFileID.close();   
  
}

void GenerateTracerMesh(long int *pktdetails) {
  
  if(myrank == 0) cout << "Generating Tracer mesh: " << flush;
  
  long int jbegin = pktdetails[1]; 
  long int jend = pktdetails[2];
  
  MyTrace_NumDrifters = jend-jbegin+1;
  
  Trace_Array = new struct Trace_Point[MyTrace_NumDrifters];
  if(Trace_Array == NULL) {
    cout << "Memory Allocation Error" << endl;
    exit(1);
  }
  
  int ij[ND];
  int counter=0;
  for(int tc=0;tc<MyTrace_NumDrifters;++tc) {
    f2ij(jbegin+tc,ij,Trace_MeshRes);
    for(int d=0;d<ND;++d)
      Trace_Array[tc].X[d] = Trace_MeshMin[d] + ij[d]*Trace_MeshDelta[d];
    Trace_Array[tc].ReleaseTime = Trace_MeshReleaseTime;
    Trace_Array[tc].Color = sqrt(Trace_Array[tc].X[0]*Trace_Array[tc].X[0]+Trace_Array[tc].X[1]*Trace_Array[tc].X[1]);//Trace_Array[tc].X[Trace_ColorDimension];
    Trace_Array[tc].LeftDomain = TestOutsideDomain(Trace_Array[tc].X, Trace_MeshReleaseTime);  
    if(!Trace_Array[tc].LeftDomain) counter++;
  }
  
  if(myrank == 0) cout << counter << " drifters in domain." << endl;
}

void ReadInTracerData(long int *pktdetails) {
  
  long int jbegin = pktdetails[1]; 
  long int jend = pktdetails[2];
  MyTrace_NumDrifters = jend-jbegin+1;
  
  ifstream Trace_InFileID(Trace_Input);
  if(!Trace_InFileID.is_open()) {
    cout << "Error opening file: \a" << Trace_Input << endl;
    exit(FAILURE);
  }
  
  char buf[LONGSTRING];
  int  numdrifters;
  Trace_InFileID.getline(buf,LONGSTRING,'\n');
  sscanf(buf, "%d\n", &numdrifters);
  if(numdrifters!=Trace_NumDrifters) {
    cout << "\aERROR: Number of drifters in " << Trace_Input << "has changed!" << endl;
    exit(1);
  }
  
  Trace_Array = new struct Trace_Point[MyTrace_NumDrifters];
  if(Trace_Array == NULL) {
    cout << "Memory Allocation Error" << endl;
    exit(1);
  }
  
  char *firstchar;
  char *eptr;
  
  for(int ii=0;ii<jbegin;++ii)
    Trace_InFileID.getline(buf,LONGSTRING,'\n');
  
  for(int ii = 0; ii < MyTrace_NumDrifters; ii++) {
    Trace_InFileID.getline(buf,LONGSTRING,'\n');
    firstchar = &buf[0];
    for(int d=0;d<ND;++d) {
      Trace_Array[ii].X[d] = strtod(firstchar,&eptr);
      firstchar = eptr;
    }
    Trace_Array[ii].ReleaseTime = strtod(firstchar,&eptr);
    firstchar = eptr;
    Trace_Array[ii].Color = strtod(firstchar,&eptr);
    Trace_Array[ii].LeftDomain = TestOutsideDomain(Trace_Array[ii].X,Trace_Array[ii].ReleaseTime); 
  }
  Trace_InFileID.close();
}

void IntegrateTracers(double t1, double t) {
  
  double tracestarttime;
  
  for(int ii = 0; ii < MyTrace_NumDrifters; ii++) {
    if(!Trace_Array[ii].LeftDomain || (Velocity_Format <= 1)) {
      if((Trace_Array[ii].ReleaseTime-t)*Time_Direction <= 0.0) {                      
        tracestarttime = Time_Direction*fmax(Time_Direction*t1, Time_Direction*Trace_Array[ii].ReleaseTime);
        if((t-tracestarttime)*Time_Direction > 0) {
          
          Integrate(Trace_Array[ii].X, tracestarttime, t);
          
          if(TestOutsideDomain(Trace_Array[ii].X,t)) 
            Trace_Array[ii].LeftDomain = TRUE; 
          
        }
      }
    }
  }  
}


void OutputTracers(double time) {    
  
 
  
  /* Write out in binary format a snapshot of the current tracers at this time for this packet */
  
  double modval;
  vector<int> activedrift;
  
  /* How many drifters do I have left in the domain at this time ?*/
  for(int ii = 0; ii < MyTrace_NumDrifters; ++ii) {
    for(int d=0;d<ND;++d) {
      if(!isfinite(Trace_Array[ii].X[d])) {
        Trace_Array[ii].LeftDomain=1;
        break;
      }
    }
    
    if(!Trace_Array[ii].LeftDomain) {
      if((Trace_Array[ii].ReleaseTime-time)*Time_Direction <= 0) {
        
        activedrift.push_back(ii);
        
      }
    }
  }
  
  int numdrift = (int)(activedrift.size());  /* That's how many */
  
   if (myrank == 0) cout << "Outputting " << numdrift << " tracers at t = " << time << endl;
  
  /* Write out coordinate of surviving tracers along each dimension */
  vector<int>::iterator iter;
  for(int d=0;d<ND;++d) {
    /* Write header */
    Trace_OutFileID.write((char *)&numdrift, sizeof(int) );
    Trace_OutFileID.write((char *)&time, sizeof(double) );
    for( iter = activedrift.begin(); iter != activedrift.end(); iter++ ) {
      if(Data_Periodic[d]) {
        modval = fmod(Trace_Array[*iter].X[d]-Data_Min[d], Data_Period[d]);
        if(modval < 0.0) {
          modval += Data_Period[d];
        }
        modval += Data_Min[d];
      }
      else {
        modval = Trace_Array[*iter].X[d];
      }
      Trace_OutFileID.write((char *)&modval, sizeof(double) );
    }
  }
  
  /* Write out color along each dimension */
  Trace_OutFileID.write((char *)&numdrift, sizeof(int) );
  Trace_OutFileID.write((char *)&time, sizeof(double) );

  for( iter = activedrift.begin(); iter != activedrift.end(); iter++ ) 
    Trace_OutFileID.write((char *)&Trace_Array[*iter].Color, sizeof(double) );
    
  
  
}

