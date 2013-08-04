// FTLE.cpp
//
// FTLE module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "ftle.h"
#include "boundary.h"
#include "parallel.h"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void OutputFTLEtoFile(int ss) {
  
  char outfile[LONGSTRING];

  sprintf(outfile,"%s%04d.raw",FTLE_Output,ss);
  
  double launchtime=Output_T1+ss*Output_TDelta;
  
  double *FTLEarray = new (nothrow) double[FTLE_BlockSize];
  if(FTLEarray == NULL) {
    cout << "Memory Allocation error: FTLEtemp" << endl;
    exit(1);
  }
  
  char infile[LONGSTRING];
  sprintf(infile, "%s%04d.ftle",FTLE_Work,ss);
  
  ifstream fin(infile,ios::binary);
  if(fin == NULL) {
    cout << "\n\aError opening " << infile << endl;
    exit(1);
  }
  
  fin.read((char *)&FTLEarray[0],FTLE_BlockSize*sizeof(double));
  fin.close();
  
  char syscommand[LONGSTRING];
  sprintf(syscommand,"rm -f %s",infile);
  system(syscommand);
  
  double *smFTLE, *Omega, *Eval;
  if(LCS_Extract) {
    /* Memory for scalar quantities for LCS extraction */
    smFTLE = new (nothrow) double[FTLE_BlockSize];
    if(smFTLE==NULL) {
      cout << "Memory Allocation error: smFTLE" << endl;
      exit(1);
    }
    Omega = new (nothrow) double[FTLE_BlockSize];
    if(Omega==NULL) {
      cout << "Memory Allocation error: Omega" << endl;
      exit(1);
    }
    Eval = new (nothrow) double[FTLE_BlockSize];
    if(Eval==NULL) {
      cout << "Memory Allocation error: Eval" << endl;
      exit(1);
    }
    
    Get_OmegaEval(FTLEarray,Omega,Eval,smFTLE);
  }
  
  int zero = 0;
  
  ofstream fout(outfile,ios::binary);
  if(fout == NULL) {
    cout << "\n\aError opening " << outfile << endl;
    exit(1);
  }
  
  fout.write((char *)&zero, sizeof(int) );
  fout.write((char *)&Time_Origin[0], 6*sizeof(int) );
  fout.write((char *)&launchtime, sizeof(double) );
  fout.write((char *)&Output_TRes, sizeof(int) );
  fout.write((char *)&Atmos_Set, sizeof(int) );
  fout.write((char *)&Atmos_Radius, sizeof(double) );
  fout.write((char *)&ss, sizeof(int) );
  fout.write((char *)&Track_Storm, sizeof(int) );
  
  /* Set minimum of FTLE grid for each dimension */
  double ftlemin[ND];
  if(Track_Storm) {
    /* Determine track coordinates by linear interpolation */
    vector<double> trackarray;
    int trackres = ReadInTrack(trackarray);
    
    double tloc;
    double tmod =  launchtime-Data_TMin;
    if(tmod > Data_TDelta * (trackres - 1)) {
      cout << "ERROR: Track Data only goes up to time " << Data_TDelta * (trackres - 1) << endl;
      exit(1);
    }
    
    int t1=(int)( tmod / Data_TDelta );
    
    if(t1==trackres-1) {
      t1=Data_TRes-2;
      tloc = 1.0;
    }
    else { 
      tloc=(tmod-t1*Data_TDelta)/Data_TDelta;
    }
    
    for(int d=0;d<ND;++d) {
      if(d<2) 
        ftlemin[d] = (1.0-tloc)*trackarray[t1*ND + d]+tloc*trackarray[(t1+1)*ND + d] - FTLE_TrackWidth[d]/2.0;
      else
        ftlemin[d] = FTLE_Min[d]; /* only x,y direction is tracked - just usual range */
    } 
  }
  else {
    for(int d=0;d<ND;++d)
      ftlemin[d] = FTLE_Min[d];
  }
  
  
  fout.write((char *)&ftlemin[0], ND * sizeof(double) );
  double ftlemax[ND];
  for(int d=0;d<ND;++d)
    ftlemax[d] = ftlemin[d] + (FTLE_Res[d]-1)*FTLE_Delta[d];
  
  fout.write((char *)&ftlemax[0], ND * sizeof(double) );
  fout.write((char *)&FTLE_Res[0], ND * sizeof(int) );
  fout.write((char *)&LCS_NumFields, sizeof(int) );
  fout.write((char *)FTLEarray,FTLE_BlockSize * sizeof(double));
  
  if(LCS_Extract) {
    fout.write((char *)Omega,FTLE_BlockSize * sizeof(double));
    fout.write((char *)Eval,FTLE_BlockSize * sizeof(double));
    fout.write((char *)smFTLE,FTLE_BlockSize * sizeof(double));
  }
  delete []FTLEarray;
  
  if(LCS_Extract) {
    delete []smFTLE;
    delete []Omega;
    delete []Eval; 
  }
  
  fout.close();
}

void CopyFTLEtoWork(long int pktnum) {

  char syscommand[LONGSTRING];
  
  sprintf(syscommand,"mv %s*_%04ld.ftle %s",FTLE_Scratch,pktnum,Path_Work);
  system(syscommand);
  

}



void CombineFTLEFiles(int ss) {
  
  char syscommand[LONGSTRING];
  
  if (myrank==0 && ss==0)
    cout << "Combining FTLE files ... " << flush;
  sprintf(syscommand,"cat %s%04d_*.ftle > %s%04d.ftle",FTLE_Work,ss,FTLE_Work,ss);
  system(syscommand);
  
  sprintf(syscommand,"rm -f %s%04d_*.ftle",FTLE_Work,ss);
  system(syscommand);
  if(myrank == 0)
    cout << "OK!" << endl;
  
}


int	ReadInTrack(vector<double>& trackarray) {
  
  char buf[LONGSTRING];
  char *eptr;
  char *firstchar;
	if(myrank == 0) cout << "Reading in Track Coordinates from " << Track_Input << " ... " << flush;
  
  /* Open storm center file */  
  ifstream Track_InFileID(Track_Input);
  if(!Track_InFileID.is_open()) {
    cout << "ERROR opening file: \a" << Track_Input << endl;
    exit(FAILURE);
  }
  
  /* Read in data */
  Track_InFileID.ignore(LONGSTRING,'\n');
  Track_InFileID.ignore(LONGSTRING,'\n');
  Track_InFileID.ignore(LONGSTRING,'\n');
  Track_InFileID.getline(buf,LONGSTRING,'\n');
  int trackres = (int)strtold(buf,NULL);
  
  trackarray.resize(trackres * ND, 0);
  
  for(int index=0;index<trackres;++index) {
    Track_InFileID.getline(buf,LONGSTRING,'\n');
    firstchar = &buf[0];
    for(int d=0;d<ND;++d) {
      trackarray[index*ND + d] = strtod(firstchar,&eptr);
      firstchar = eptr;
    }
  }
  Track_InFileID.close();
  
  if(myrank == 0) cout << "Read in " << trackres << " points." << endl;
  
  
  return(trackres);
}


void ReadInFTLESlide(int ss, long int packetID) {
  
  /* Function reads in FTLE data structures for slide specified by ss into FTLE_Array */
  char ftlebinfile[LONGSTRING];
  sprintf(ftlebinfile, "%s%04d_%04ld.man",FTLE_Scratch,FTLE[ss].slidenumber,packetID);
  ifstream FTLE_BinFileID(ftlebinfile,ios::binary);
  if(!FTLE_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << ftlebinfile << endl;
    exit(FAILURE);
  }
  
  /* Read FTLE_Array information from file */
  
  FTLE_BinFileID.read((char *)&FTLE_Pts[0], sizeof(struct FTLE_Point) * FTLE[ss].blocksize);
  
  int numdrifters = 0;
  FTLE_BinFileID.read((char *)&numdrifters, sizeof(int));
  int *itmp = new int[numdrifters]; 
  FTLE_BinFileID.read((char *)&itmp[0], numdrifters*sizeof(int));
  
  double *dtmp = new double[ND*numdrifters]; 
  FTLE_BinFileID.read((char *)&dtmp[0], ND*numdrifters*sizeof(double));
  
  FTLE_Dfts.clear();
  for(int ii=0;ii<numdrifters;++ii) {
    int key = itmp[ii];
    for(int d=0;d<ND;++d)
      FTLE_Dfts[key].X[d]=dtmp[ii*ND+d];
  }
  delete []itmp;
  delete []dtmp;
  
  FTLE_BinFileID.close();
  
}

void WriteOutFTLESlide(int ss, long int packetID) {
  /* Writes contents of FTLE_Array to binary file for slide ss */
  
  char ftlebinfile[LONGSTRING];
  sprintf(ftlebinfile, "%s%04d_%04ld.man",FTLE_Scratch,FTLE[ss].slidenumber,packetID);
  
  ofstream FTLE_BinFileID(ftlebinfile,ios::binary);
  if(!FTLE_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << ftlebinfile << endl;
    exit(FAILURE);
  }
  
  FTLE_BinFileID.write((char *)&FTLE_Pts[0], sizeof(FTLE_Pts[0]) * FTLE[ss].blocksize);
  
  int numdrifters = FTLE_Dfts.size();
  FTLE_BinFileID.write((char *)&numdrifters, sizeof(int));
  
  map<int,struct FTLE_Drift>::iterator it;
  for(it = FTLE_Dfts.begin();it != FTLE_Dfts.end() ; ++it)
    FTLE_BinFileID.write((char *)&(*it).first, sizeof(int));
  for(it = FTLE_Dfts.begin();it != FTLE_Dfts.end() ; ++it)
    FTLE_BinFileID.write((char *)&(*it).second, ND*sizeof(double));
  FTLE_BinFileID.close();
  
}


void AllocateMemoryForFTLE(long int *pktdetails) {
  
  long int j1,j2;
  long int packetID = pktdetails[0];
  long int jbegin = pktdetails[1]; 
  long int jend = pktdetails[2];
  
  long int ssbegin = (long int)(jbegin/FTLE_BlockSize);
  long int ssend = (long int)(jend/FTLE_BlockSize);
  
  
  FTLE_Pts = new (nothrow) struct FTLE_Point[FTLE_BlockSize];
  if(FTLE_Pts==NULL) {
    cout << "Memory Allocation error" << endl;
    exit(1);
  }
  
  FTLE = new (nothrow) struct FTLE_Slide[MyOutput_TRes];
  if(FTLE==NULL) {
    cout << "Memory Allocation error" << endl;
    exit(1);
  }
  vector<double> trackarray;
  int trackres = 0;
  if(Track_Storm) {
    /* Determine track coordinates by linear interpolation */
    trackres = ReadInTrack(trackarray);
  }
  
  /* Initiate slide information for each slide */
  for(int ss = 0; ss < MyOutput_TRes; ss++) {
    
    ss==0 ? (j1 = jbegin-(ssbegin*FTLE_BlockSize)) : (j1=0);
    ss==(MyOutput_TRes-1) ? (j2 = jend-(ssend*FTLE_BlockSize)) : (j2=FTLE_BlockSize-1);
    
    FTLE[ss].blocksize = j2-j1+1;
    FTLE[ss].numdrifters = 0;
    FTLE[ss].slidenumber = ssbegin+ss;
    FTLE[ss].completestatus = 2;
    FTLE[ss].launchtime = Output_T1 + (ssbegin+ss) * Output_TDelta;
    
    FTLE[ss].stoptime = FTLE[ss].launchtime + Time_Direction * FTLE_IntTLength;
    if(!Data_TPeriodic && (Velocity_Format > 1)) {
      FTLE[ss].stoptime = fmin(FTLE[ss].stoptime, Data_TMax);
      FTLE[ss].stoptime = fmax(FTLE[ss].stoptime, Data_TMin);
    }
    
    
    /* Set minimum of FTLE grid for each dimension */
    
    if(Track_Storm) {
      double tloc;
      double tmod =  FTLE[ss].launchtime-Data_TMin;
      if(tmod > Data_TDelta * (trackres - 1)) {
        cout << "ERROR: Track Data only goes up to time " << Data_TDelta * (trackres - 1) << endl;
        exit(1);
      }
      
      int t1=(int)( tmod / Data_TDelta );
      
      if(t1==trackres-1) {
        t1=Data_TRes-2;
        tloc = 1.0;
      }
      else { 
        tloc=(tmod-t1*Data_TDelta)/Data_TDelta;
      }
      
      for(int d=0;d<ND;++d) {
        
        if(d<2) 
          FTLE[ss].ftlemin[d] = (1.0-tloc)*trackarray[t1*ND + d]+tloc*trackarray[(t1+1)*ND + d] - FTLE_TrackWidth[d]/2.0;
        else
          FTLE[ss].ftlemin[d] = FTLE_Min[d]; /* only x,y direction is tracked - just usual range */
      } 
    }
    else {
      for(int d=0;d<ND;++d)
        FTLE[ss].ftlemin[d] = FTLE_Min[d];
    }
    
    
    /*Set information about FTLE points and drifters */
    for(int f = 0; f < FTLE[ss].blocksize; ++f) {
      
      FTLE_Pts[f].HaveFTLE = 0;
      FTLE_Pts[f].FTLEValue = 0.0;
      
      int ij0[ND];
      f2ij(j1+f,ij0,FTLE_Res);
      /*test if nbrs are outside domain*/
      for(int d1=0;d1<ND;++d1) {
        
        double nbr0[ND];
        double nbr1[ND];
        
        /* Generate two neighbors, one dimension at a time */
        for(int d=0;d<ND;++d) {
          nbr0[d]=FTLE[ss].ftlemin[d]+ij0[d]*FTLE_Delta[d];
          nbr1[d]=nbr0[d];
        }
        nbr0[d1]-=FTLE_Delta[d1]; 
        nbr1[d1]+=FTLE_Delta[d1];
        
        
        /* Test if neighbors are outside domain */
        /* If outside, we cannot compute FTLE; set FTLE to 0 */
        /* If not outside, make a drifter */
        if( TestOutsideDomain(nbr0,FTLE[ss].launchtime) || TestOutsideDomain(nbr1,FTLE[ss].launchtime) ) {
          FTLE_Pts[f].HaveFTLE = 1;
          FTLE_Pts[f].FTLEValue = 0.0;
          continue;
        }
      } 
      
      /* Get absolute indices of neighbors, if do not exist, then we must create them. */
      
      if(!FTLE_Pts[f].HaveFTLE) {
        int fn;
        int ijn[ND];
        
        for(int d1=0;d1<ND;++d1) {
          for(int d=0;d<ND;++d)
            ijn[d]=ij0[d]+1;  /*Dft grid is one offset from pts grid */
          
          /*left neighbor */
          ijn[d1]=ij0[d1];
          
          fn=ij2f(ijn,FTLE_DftRes);
          for(int d=0;d<ND;++d)
            FTLE_Dfts[fn].X[d]=FTLE[ss].ftlemin[d]+(ijn[d]-1)*FTLE_Delta[d];
          
          /*right neighbor */
          ijn[d1]=ij0[d1]+2;  /*Dft grid is one offset from pts grid */
          
          fn=ij2f(ijn,FTLE_DftRes);
          for(int d=0;d<ND;++d)
            FTLE_Dfts[fn].X[d]=FTLE[ss].ftlemin[d]+(ijn[d]-1)*FTLE_Delta[d];
          
        }
      }
    } 
    WriteOutFTLESlide(ss,packetID);
  }  
}

void FreeMemoryForFTLE() {
  
  delete []FTLE_Pts;  
  delete []FTLE;
  
}

void	GetFTLE(int fn, double FTLEtime, int j1, int j2) { 
  /* Compute the FTLE value for all points for which fn is a neighbor */
  /* Remember that fn is on the drifter grid, and the FTLE points are on the smaller points grid */
  
  // cout << myrank  <<  ": GetFTLE fn = " << fn << endl;
  
  int ij[ND];
  int ijn[ND];
  int f;
  
  f2ij(fn,ijn,FTLE_DftRes);
  
  for(int d=0;d<ND;++d) {
    
    for(int di=0;di<ND;++di) 
      ij[di]=ijn[di]-1;  /* ij is in pts grid coordinates */
    
    ij[d]=ijn[d]-2;  /* left neighbor */
    
    int skipflag=0;
    for(int di=0;di<ND;++di) {
      if( (ij[di]<0) || (ij[di] > (FTLE_Res[di]-1)) ) {
        skipflag=1;
        break;
      }
    }
    
    if(!skipflag) {
      
      f = ij2f(ij,FTLE_Res); 
      
      if(f>=j1 && f<=j2) {
        if(!(FTLE_Pts[f-j1].HaveFTLE)) {
          FTLE_Pts[f-j1].FTLEValue=GetFTLEForPoint(f,FTLEtime);
          FTLE_Pts[f-j1].HaveFTLE=TRUE;	
        }
      }
    }
    
    ij[d]=ijn[d]; /* right neighbor */
    
    skipflag=0;
    for(int di=0;di<ND;++di) {
      if( (ij[di]<0) || (ij[di] > (FTLE_Res[di]-1)) ) {
        skipflag=1;
        break;
      }
    }
    
    if(!skipflag) {
      
      f = ij2f(ij,FTLE_Res);  
      
      if(f>=j1 && f<=j2) {
        if(!(FTLE_Pts[f-j1].HaveFTLE)) {
          FTLE_Pts[f-j1].FTLEValue=GetFTLEForPoint(f,FTLEtime);
          FTLE_Pts[f-j1].HaveFTLE=TRUE;	
        }
      }
    }
  }
}

void Filter(double FTLEarray[], double smFTLE[], double sigma) {
  
  int ij[ND];
  
  /* Build filter kernel */
  // flength is odd by construction.
  int midpt = (int)floor(sigma*3.0);
  int flength = (2 * midpt) + 1;
  double *fvals = new (nothrow) double[flength];
  
  double fsum = 0.0;
  for(int w=0;w<flength;++w) {
    fvals[w] = exp(-(w-midpt)*(w-midpt)/(2*sigma*sigma));
    fsum+=fvals[w];
  }
  for(int w=0;w<flength;++w) 
    fvals[w]/=fsum;
  
  /* Perform Sequence of 1D Convolutions */
  int ijnb[ND];
  for(int f=0;f<FTLE_BlockSize; ++f) {
    double newval = 0.0;
    fsum = 0.0;
    f2ij(f,ij,FTLE_Res);
    for(int d=0;d<ND;++d) {
      for(int w=0; w<flength; ++w) {
        int p=ij[d]-midpt+w;
        if( (p > -1) && (p < FTLE_Res[d]) ) {
          for(int n=0;n<ND;++n)
            ijnb[n]=ij[n];
          ijnb[d] = ij[d]-midpt+w;
          newval += (fvals[w]*FTLEarray[ij2f(ijnb,FTLE_Res)]);
          fsum += fvals[w];
        }
      }
    }
    smFTLE[f] = newval/fsum;
  }
  
  delete []fvals;
  
}

void copy_vi(int *v1, int *v2, int n) {
  
  for(int d=0;d<n;++d)
    v2[d] = v1[d];
  
  return;
  
}

int fn(int *c, int *p, int *m, int d) {
  
  copy_vi(c,p,ND);
  copy_vi(c,m,ND);
  if(p[d] < FTLE_Res[d]-1) p[d]+=1;
  if(m[d] > 0) m[d]-=1;
  
  if(c[d] == 0 || c[d] == FTLE_Res[d]-1 )
    return 1;
  else
    return 2;
}

double Get_D(double *field, int f, int d1, int d2) {
  
  // diagonal elements
  
  int ijc[ND];
  int ijp[ND];
  int ijm[ND];
  
  f2ij(f,ijc,FTLE_Res);
  
  // shift if on the boundary
  for(int d=0;d<ND;++d) {
    if(ijc[d] == 0) 
      ijc[d]+=1;
    else if(ijc[d] == FTLE_Res[d]-1) 
      ijc[d]-=1;
  }
  
  // diagonal terms
  if(d1==d2) {
    
    copy_vi(ijc,ijm,ND);
    ijm[d1]-=1;
    
    copy_vi(ijc,ijp,ND);
    ijp[d1]+=1;
    
    double cval = field[ij2f(ijc,FTLE_Res)];
    double pval = field[ij2f(ijp,FTLE_Res)];
    double mval = field[ij2f(ijm,FTLE_Res)];
    return ( pval - 2.0*cval + mval) / (FTLE_Delta[d1] * FTLE_Delta[d1]);
  }
  else { //off-diagonal terms
    
    copy_vi(ijc,ijp,ND);
    copy_vi(ijc,ijm,ND);
    
    ijp[d1]+=1;
    ijp[d2]+=1;
    ijm[d1]-=1;
    ijm[d2]+=1;
    double dp1 = ( field[ij2f(ijp,FTLE_Res)] - field[ij2f(ijm,FTLE_Res)] ) / (2.0 * FTLE_Delta[d1]);
    
    ijp[d2]-=2;
    ijm[d2]-=2;
    double dm1 = ( field[ij2f(ijp,FTLE_Res)] - field[ij2f(ijm,FTLE_Res)] ) / (2.0 * FTLE_Delta[d1]);
    
    return ( dp1 - dm1 ) / ( 2.0 * FTLE_Delta[d2] );
  }
  
}


double Get_D(double *field, int fc, int d) {
  
  int ijc[ND];
  int ijp[ND];
  int ijm[ND];
  
  f2ij(fc,ijc,FTLE_Res);
  
  int numspaces = fn(ijc,ijp,ijm,d);
  
  double pval = field[ij2f(ijp,FTLE_Res)];
  double mval = field[ij2f(ijm,FTLE_Res)];
  
  return ( pval - mval ) / (numspaces * FTLE_Delta[d]);
}

int Get_OmegaEval(double FTLEarray[], double omega[], double mineval[], double smFTLE[]) {
  
  
  if(myrank == 0)
    cout << "Extracting LCS " << flush; 
  
  /* Prepare the FTLE field for derivatives by smoothing */
  
  if(Filter_Width) {
    Filter(FTLEarray,smFTLE,Filter_Width);
  }
  else {
    for(int f=0;f<FTLE_BlockSize;++f) 
      smFTLE[f] = FTLEarray[f];
  }
  
  /* Compute Gradient, Hessian, Eigenvalues, Eigenvectors, and Dot product of mineig and grad. */
  
  double grad[ND];
  double ddval;
  int statusstep = FTLE_BlockSize/10;
  
  gsl_matrix *ddx = gsl_matrix_alloc(ND,ND);
  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(ND);
  gsl_vector *eval = gsl_vector_alloc (ND);
  gsl_matrix *gslevec = gsl_matrix_alloc(ND,ND);
  
  for(int f=0;f<FTLE_BlockSize;++f) {
    if(myrank == 0) {
      if(!(f % statusstep)) cout << "." << flush;
    }
    
    for(int d=0;d<ND;++d) {
      
      //get gradient
      
      grad[d] = Get_D(smFTLE,f,d);
      
      //get Hessian
      for(int ds=0;ds<=d;++ds) {
        
        ddval = Get_D(smFTLE, f, d, ds);
        
        if(!isfinite(ddval)) {
          cout << "Error: value in Hessian is " << ddval << endl;
        } 
        
        gsl_matrix_set(ddx,d,ds,ddval);
        
      }
    }
    
    // compute eigenvalues 
    
    int status = gsl_eigen_symmv(ddx, eval, gslevec, work);
    
    if( status != GSL_SUCCESS ) gsl_strerror(status);
    
    int minindex = gsl_vector_min_index(eval);
    
    mineval[f] = gsl_vector_get(eval,minindex);
    
    
    double evec[ND];
    for(int d=0;d<ND;++d)
      evec[d] = gsl_matrix_get(gslevec,d , minindex);
    
    
    if( evec[1] < 0.0 )
      for(int d=0;d<ND;++d)
        evec[d] *= (-1.0);
    
    double dot = 0.0;
    for(int d=0;d<ND;++d)
      dot += (evec[d]*grad[d]);
    omega[f]=dot;
    
  }
  
  gsl_matrix_free(ddx);
  gsl_matrix_free(gslevec);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(work);
  
  
  return 0;
}




void ComputeFTLE(int ss, long int j1, long int packetID) {
  /* Compute and Write FTLE values for slide ss to a binary file */
  
  for(int f=0 ; f < FTLE[ss].blocksize; f++) {
    
    if(!(FTLE_Pts[f].HaveFTLE)) 
      FTLE_Pts[f].FTLEValue = GetFTLEForPoint(j1+f, FTLE_IntTLength);
    
  }
  
  char outfile[LONGSTRING];
  sprintf(outfile, "%s%04d_%04ld.ftle",FTLE_Scratch,FTLE[ss].slidenumber,packetID);
  
  ofstream fout(outfile,ios::binary);
  if(fout == NULL) {
    cout << "\n\aError opening " << outfile << endl;
    exit(1);
  }
  
  double *tmpd = new double[FTLE[ss].blocksize];
  for(int ii=0;ii<FTLE[ss].blocksize;++ii) 
    tmpd[ii]=FTLE_Pts[ii].FTLEValue;
  
  fout.write((char *)&tmpd[0],FTLE[ss].blocksize*sizeof(double));
  
  fout.close();
  
  /* Delete the old slide file */
  char syscommand[LONGSTRING];
  sprintf(syscommand,"rm -f %s%04d_%04ld.man",FTLE_Scratch,FTLE[ss].slidenumber,packetID);
  system(syscommand);
  
}

double Hamiltonian(double X[ND]) {

  double eterm = (-1.0+exp(-7.0*(0.7-cos(X[0]))));
  return ( 0.5*X[1]*X[1] + 1.0/1400.0*eterm*eterm );

}

double	GetFTLEForPoint(int f0, double FTLEtime) {
  
  if(FTLEtime < 0.001 * FTLE_IntTLength) {
    return 0.0;
  }
  else {
    double pts[ND][2][ND];
    double lambda;
    
    BuildDifferenceMatrix(f0, pts);
    
    
    if(FTLE_Type == 0) {
      
      //Define matrix Dphi
      double Dphi[ND][ND];
      
      for(int d1=0;d1<ND;++d1) {
        for(int d2=0;d2<ND;++d2) {
          Dphi[d1][d2] = (pts[d2][1][d1] - pts[d2][0][d1]) / (2 * FTLE_Delta[d2]);
        }
      }
      
      // Delta = Dphi^T * Dphi
      double Delta[ND][ND];
      for(int ii = 0; ii < ND; ii++) {
        for(int jj = 0; jj <= ii; jj++) {
          Delta[ii][jj] = 0;
          for(int mm = 0; mm < ND; mm++)
            Delta[ii][jj] += Dphi[mm][ii] * Dphi[mm][jj];
        }
      }
      
      /* GSL only needs lower half of symmetric matrix */
      gsl_matrix *Delta_gsl = gsl_matrix_alloc(ND,ND);
      for(int ii = 0; ii < ND; ii++) {
        for(int jj = 0; jj <= ii; jj++) {
          gsl_matrix_set(Delta_gsl,ii,jj,Delta[ii][jj]);
        }
      }
      
      // Get eigenvalues of Delta
      gsl_eigen_symm_workspace *work = gsl_eigen_symm_alloc(ND);
      gsl_vector *eval = gsl_vector_alloc (ND);
      int status = gsl_eigen_symm (Delta_gsl, eval, work);
      
      if( status != GSL_SUCCESS ) gsl_strerror(status);
      
      lambda =  gsl_vector_max(eval);
      
      //Clean up
      gsl_matrix_free(Delta_gsl);
      gsl_vector_free(eval);
      gsl_eigen_symm_free(work);
      
      if(lambda>1.0) {  /* lambda > 0.0   ?????????*/
        return 0.5*log(lambda)/FTLEtime;   /* Option to divide by FTLE_IntTLength or FTLEtime here */
        /*return 0.5*log(lambda);*/
      }
      else {
        return 0.0;
      }
    }
    else if(FTLE_Type == 1) {
    
#if ND==2
      const double saddle = 0.000714276;
      
      for(int d1=0;d1<ND;++d1) {
        int out[2];
        for(int nb=0;nb<2;++nb) {
          if( (pts[d1][nb][0]< 0.0) || (pts[d1][nb][0]>M_PI) || (Hamiltonian(pts[d1][nb]) > saddle) )
            out[nb] = 1;
          else
            out[nb] = 0;
        }
        if(out[0]-out[1])
          return 1.0;
      }
      
      return 0.0;
      
#elif ND==3

    const double saddle = 0.55;
      
      for(int d1=0;d1<ND;++d1) {
        int out[2];
        for(int nb=0;nb<2;++nb) {
          if( (pts[d1][nb][0]< -1.0) || (pts[d1][nb][0]>saddle)  )
            out[nb] = 1;
          else
            out[nb] = 0;
        }
        if(out[0]-out[1])
          return 1.0;
      }
      
      return 0.0;
      
#elif ND==4

      for(int d1=0;d1<ND;++d1) {
        int out[2];
        for(int nb=0;nb<2;++nb) {
          double xval = pts[d1][nb][0];
          double yval = pts[d1][nb][1];
          double distancefromsun = (xval+FBP_mu)*(xval+FBP_mu)+yval*yval;
          
          if( distancefromsun > 1 )
            out[nb] = 1;
          else
            out[nb] = 0;
        }
        if(out[0]-out[1])
          return 1.0;
      }
      
      return 0.0;

#endif

      
    }
  }
  return 0.0;
}






