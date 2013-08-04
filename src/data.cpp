// data.cpp
//
// data module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include <cmath>
#include <vector>
using namespace std;
#include "data.h"

/* Add new modules for converting ASCII to binary data here.
Modules must obey the following standard:

1) ASCII data contained in tecfilename where
sprintf(tecfilename,"%s.dat",Data_Input);

2)  Velocity data for each time frame is written to a separate binary file named binfilename where
sprintf(binfilename,"%sf%04d.new",Data_Scratch,tt);

3) Binary is output in a single fwrite 
fwrite( U , sizeof(double) , ND * Data_BlockSize, Data_BinFileID);
where U[Data_BlockSize*ND] is accessed by U[ijd(i,j,k,d)];  
*/

void ReadInDataNestParameters() {
  
  int FirstFrameMod;
  char binfilename[LONGSTRING];
  
  FirstFrameMod = Data_FirstFrame % Data_TRes;
  if(FirstFrameMod < 0) FirstFrameMod += Data_TRes;
  
  sprintf(binfilename,"%sf%04d_0000.new",Data_Scratch,FirstFrameMod);
  
  ifstream Data_BinFileID(binfilename);
  if(!Data_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << binfilename << endl;
    exit(FAILURE);
  }
  
  Data_BinFileID.read((char *)&Data_Min[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_Max[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_Res[0], ND * sizeof(int));
  Data_BinFileID.read((char *)&Data_Delta[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_BlockSize, sizeof(int));
  Data_BinFileID.close();
  
  
}


void ReadInDataParameters() {
  
  int FirstFrameMod;
  char binfilename[LONGSTRING];
  
  FirstFrameMod = Data_FirstFrame % Data_TRes;
  if(FirstFrameMod < 0) FirstFrameMod += Data_TRes;
  
  sprintf(binfilename,"%sf%04d.new",Data_Scratch,FirstFrameMod);
  
  ifstream Data_BinFileID(binfilename);
  if(!Data_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << binfilename << endl;
    exit(FAILURE);
  }
  
  Data_BinFileID.read((char *)&Data_Min[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_Max[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_Res[0], ND * sizeof(int));
  Data_BinFileID.read((char *)&Data_Delta[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_Period[0], ND * sizeof(double));
  Data_BinFileID.read((char *)&Data_BlockSize, sizeof(int));
  Data_BinFileID.close();
  cout << "blocksize = " << Data_BlockSize << endl;
  
  Data_HeaderSize = 4*ND*sizeof(double) + (ND+1)*sizeof(int);
  
}


void ReadInNonUniformGrids() {
  
  char buf[LONGSTRING];
  char gridfilename[LONGSTRING];
  
  for(int d=0;d<ND;++d) {
    if(Data_NonUniformGrid[d]) {
      if(myrank == 0) cout << "Reading nonuniform grid " << d << endl;
      
      
      if(strcmp(Path_Input,"pwd")) {
        sprintf(gridfilename,"%s%sgrid%04d.dat",Path_Input,Data_InFile,d);
      }
      else {    
        sprintf(gridfilename,"%sgrid%04d.dat",Data_InFile,d);
      }
      
      ifstream Grid_FileID;
      Grid_FileID.open(gridfilename);
      if(!Grid_FileID.is_open()) {
        cout << "\n Error opening file: \a" << gridfilename << endl;
        exit(FAILURE);
      }
      
      Grid_FileID.getline(buf,LONGSTRING,'\n');
      int numpts = (int)strtod(buf,NULL);
      
      if(myrank == 0) cout << "grid " << d << " file has " << numpts << " points." << endl;
      
      Data_Grid[d] = new (nothrow) double[numpts];
      if(Data_Grid[d] == NULL) {
        cout << "ERROR: Cannot create memory for Data_Grid" << endl;
        exit(1);
      }
      
      for(int pt=0;pt<numpts;++pt) {
        Grid_FileID.getline(buf,LONGSTRING,'\n');
        Data_Grid[d][pt] = strtod(buf,NULL);
      }
      
      Grid_FileID.close();
      
    }
  }
}

void ReadInASCIIVelocityData(int format, int numprocs) {
  
  switch ( format ) {
    case 1 : 
      if(myrank==0) ReadInASCIITecplotVelocityData();
      break;
    case 2 : 
      for(int nf=myrank;nf<Data_NumInputFiles; nf+=numprocs) {
        ReadInASCIIMultiLeanVelocityData(nf);
      }
      break;
    case 3 : 
      for(int nf=myrank;nf<Data_NumInputFiles; nf+=numprocs) {
        ReadInASCIIMultiLeanVelocityNestData(nf,Nest_NumNests);
      }
      break;
    default : 
      if(myrank==0) ReadInASCIITecplotVelocityData();
  }
  
}


void ReadInASCIIMultiLeanVelocityNestData(int nf, int Nest_NumNests) {
  
  /*Lean ASCII files contain a general header line and three comment lines for each dimension 
  They also contain a ZONE heading for each time frame.*/
  
  char tecfilename[LONGSTRING];
  char binfilename[LONGSTRING];
  char buf[LONGSTRING];
  char *eptr;
  char *firstchar;
  
  sprintf(buf,"%d : Converting data frame %d to binary.\n",myrank,nf);
  cout << buf << flush;
  
  for(int nn=0;nn<Nest_NumNests;++nn) {  
    
    int nestID = Nest_List[nn];
    
    /* Open ASCII file for reading */
    sprintf(tecfilename,"%s%04d_%04d.dat",Data_Input,nf,nestID);
    
    ifstream Data_InFileID;
    Data_InFileID.open(tecfilename);
    if(!Data_InFileID.is_open()) {
      cout << "\n Error opening file: \a" << tecfilename << endl;
      exit(FAILURE);
    }
    
    /* Strip off Title */
    Data_InFileID.getline(buf, LONGSTRING, '\n');        
    
    double dmin[ND];
    double dmax[ND];
    int dres[ND];
    double ddelta[ND];
    
    /* Get Info about Data */
    for(int d = 0; d < ND; ++d) {
      Data_InFileID.getline(buf,LONGSTRING,'\n');
      dmin[d] = strtod(&buf[13],&eptr);
      Data_InFileID.getline(buf,LONGSTRING,'\n');
      dmax[d] = strtod(&buf[13],&eptr);
      Data_InFileID.getline(buf,LONGSTRING,'\n');
      dres[d] = (int)strtod(&buf[13],&eptr);
      ddelta[d] = (dmax[d]-dmin[d])/(dres[d]-1);
    }
    
    /* Strip off Zone Header */
    Data_InFileID.getline(buf, LONGSTRING, '\n');          
    
    if( strncmp(buf,"Z",1) ) {
      cout << "ERROR: Could not find zone heading in " << tecfilename << endl;
      exit(1);
    }
    /* Read in velocity data */
    
    int dblocksize = 1;
    for(int d=0;d<ND;++d)
      dblocksize *= dres[d];
    
    vector<double> U(dblocksize*ND,0.0); 
    for(int index=0;index<dblocksize;++index) {
      Data_InFileID.getline(buf,LONGSTRING,'\n');
      firstchar = &buf[0];
      for(int d=0;d<ND;++d) {
        U[d*dblocksize+index] = strtod(firstchar,&eptr);
        firstchar = eptr;
      }
    }
    /* Close Input File */
    Data_InFileID.close();
    
    /* Open binary file for writing */
    sprintf(binfilename,"%sf%04d_%04d.new",Data_Work,nf,nestID);
    
    ofstream Data_BinFileID;
    Data_BinFileID.open(binfilename,ios::out|ios::trunc|ios::binary);
    if(!Data_BinFileID.is_open()) {
      cout << "\n Error opening file: \a" << binfilename << endl;
      exit(FAILURE);
    }
    
    /* Write out velocity data */
    Data_BinFileID.write( (char *)&dmin[0] , sizeof(double)*ND);
    Data_BinFileID.write( (char *)&dmax[0] , sizeof(double)*ND);
    Data_BinFileID.write( (char *)&dres[0] , sizeof(int)*ND);
    Data_BinFileID.write( (char *)&ddelta[0] , sizeof(double)*ND);
    Data_BinFileID.write( (char *)&dblocksize , sizeof(int));
    Data_BinFileID.write( (char *)&U[0] , sizeof(double)*ND*dblocksize);
    
    /* Close Output file */
    Data_BinFileID.close();
    
  }
}

void ReadInASCIIMultiLeanVelocityData(int nf) {
  
  /*Lean ASCII files contain a general header line and three comment lines for each dimension 
  They also contain a ZONE heading for each time frame.*/
  
  char tecfilename[LONGSTRING];
  char binfilename[LONGSTRING];
  char buf[LONGSTRING];
  char *eptr;
  char *firstchar;
  
  sprintf(buf,"%d : Converting data frame %d to binary.\n",myrank,nf);
  cout << buf << flush;
  
  /* Open ASCII file for reading */
  sprintf(tecfilename,"%s%04d.dat",Data_Input,nf);
  
  ifstream Data_InFileID;
  Data_InFileID.open(tecfilename);
  if(!Data_InFileID.is_open()) {
    cout << "\n Error opening file: \a" << tecfilename << endl;
    exit(FAILURE);
  }
  
  /* Strip off TITLE */
  Data_InFileID.getline(buf,LONGSTRING,'\n');
  
  double dmin[ND];
  double dmax[ND];
  int dres[ND];
  double ddelta[ND];
  double dperiod[ND];
  
  for(int d=0;d<ND;++d) {
    Data_InFileID.ignore(256,'=');
    Data_InFileID >> dmin[d];
    
    Data_InFileID.ignore(256,'=');
    Data_InFileID >> dmax[d];
    
    Data_InFileID.ignore(256,'=');
    Data_InFileID >> dres[d];
    cout << dmin[d] << "  " << dmax[d] << "   " << dres[d] << endl;
  }
  
  
  
  Data_InFileID.getline(buf,LONGSTRING,'\n');
  
  int dblocksize = 1;
  for(int d=0;d<ND;++d)
    dblocksize *= dres[d];
  
  for(int d=0;d<ND;++d) {
    
    if(dmin[d] >= dmax[d]) {
      if(myrank == 0) cout << "Parameter Error: Data_Min" << d << " >= " <<  "Data_Max" << d << " = " << dmin[d] << endl;
      exit(1);
    }
    if(dres[d] < 2) {
      if(myrank == 0) cout << "Parameter Error: Data_Res%d < 2\a" << endl;
      exit(1);
    }
    
    if(!Data_NonUniformGrid[d]) {
      ddelta[d] = (dmax[d]-dmin[d])/(dres[d] - 1);
    }
    if(Data_Periodic[d]) 
      dperiod[d] = dres[d]*(dmax[d] - dmin[d])/(dres[d]-1);
    else 
      dperiod[d] = 0.0;
  }
  
  /* Strip off Zone Header */
  Data_InFileID.getline(buf, LONGSTRING, '\n');          
  
  if( strncmp(buf,"Z",1) ) {
    cout << "ERROR: Could not find zone heading in " << tecfilename << endl;
    exit(1);
  }
  /* Read in velocity data */
  
  /* Assign temporary memory for U*/
  double *U = new (nothrow) double[ND * dblocksize];
  if(U==NULL) {
    cout << "Memory Allocation error" << endl;
    exit(1);
  }
  
  for(int index=0;index<dblocksize;++index) {
    Data_InFileID.getline(buf,LONGSTRING,'\n');

    firstchar = &buf[0];
    for(int d=0;d<ND;++d) {
      U[d*dblocksize+index] = strtod(firstchar,&eptr);
      firstchar = eptr;
      //cout << U[d*dblocksize+index] << "  " << flush;
    }
    //cout << endl;
  }
  /* Close Input File */
  Data_InFileID.close();
  
  /* Open binary file for writing */
  sprintf(binfilename,"%sf%04d.new",Data_Work,nf);
  
  ofstream Data_BinFileID;
  Data_BinFileID.open(binfilename,ios::out|ios::trunc|ios::binary);
  if(!Data_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << binfilename << endl;
    exit(FAILURE);
  }
  
  /* Write out velocity data */
  Data_BinFileID.write( (char *)&dmin[0] , sizeof(double)*ND);
  Data_BinFileID.write( (char *)&dmax[0] , sizeof(double)*ND);
  Data_BinFileID.write( (char *)&dres[0] , sizeof(int)*ND);
  Data_BinFileID.write( (char *)&ddelta[0] , sizeof(double)*ND);
  Data_BinFileID.write( (char *)&dperiod[0] , sizeof(double)*ND);
  Data_BinFileID.write( (char *)&dblocksize , sizeof(int));
  Data_BinFileID.write( (char *)&U[0] , sizeof(double)*ND*dblocksize);
  
  /* Write out velocity data */
  Data_BinFileID.write( (char *)U , sizeof(double)*ND*dblocksize);
  
  /* Close Output file */
  Data_BinFileID.close();
  
  delete []U;
}

void ReadInASCIITecplotVelocityData() {
  
  int DataFile_Headers = 1;
  char tecfilename[LONGSTRING];
  char binfilename[LONGSTRING];
  char buf[LONGSTRING];
  char *eptr;
  char *firstchar;
  
  
  /* Assign temporary memory for U*/
  double *U = new (nothrow) double[ND * Data_BlockSize];
  if(U==NULL) {
    cout << "Memory Allocation error" << endl;
    exit(1);
  }
  
  /* Open ASCII file for reading */
  sprintf(tecfilename,"%s.dat",Data_Input);
  
  ifstream Data_InFileID(tecfilename);
  if(!Data_InFileID) {
    cout << "\n Error opening file: \a" << tecfilename << endl;
    exit(FAILURE);
  }
  
  cout << "Reading in ASCII Tecplot Velocity Data from " << tecfilename << " ..." << endl;
  
  /* Read in ASCII data and write to binary file one slide of data at a time*/
  
  /* Strip off TITLE */
  for(int i = 0; i < DataFile_Headers; ++i) {
    Data_InFileID.getline(buf, LONGSTRING, '\n');
    cout << buf << endl;
  }
  
  for(int tt = 0; tt < Data_TRes; tt++) {
    /* Strip off Zone Header */
    cout << "tt = " << tt << endl;
    Data_InFileID.getline(buf, LONGSTRING, '\n');       
    cout << buf << endl;
    if( strncmp(buf,"Z",1) ) {
      cout << "ERROR: Data does not have correct resolution in " << tecfilename << endl;
      exit(1);
    }
    
    /* Read in velocity data */
    cout << "dbsz = " << Data_BlockSize << endl;
    for(int index=0;index<Data_BlockSize;++index) {
      cout << index << endl;
      Data_InFileID.getline(buf,LONGSTRING,'\n');
      firstchar = &buf[0];
      for(int d=0;d<ND;++d) {
        cout << strtod(firstchar,&eptr) << " ";
        firstchar = eptr;
      }
      for(int d=0;d<ND;++d) {
        U[d*Data_BlockSize + index] = strtod(firstchar,&eptr);
        cout << U[d*Data_BlockSize + index];
        firstchar = eptr;
      }
      cout << endl;
    }
    
    /* Open binary file for writing */
    sprintf(binfilename,"%sf%04d.new",Data_Work,tt);
    cout << "Writing zone to " << binfilename << endl;
    
    ofstream Data_BinFileID(binfilename);
    if(!Data_BinFileID.is_open()) {
      cout << "\n Error opening file: \a" << binfilename << endl;
      exit(FAILURE);
    }    
    /* Write out velocity data */
    Data_BinFileID.write( (char *)U ,ND * sizeof(double)*Data_BlockSize);
    
    /* Close Output file */
    Data_BinFileID.close();
  }
  
  /* Close Input File */
  Data_InFileID.close();
  
  delete []U;
}

void LoadFirstDataFrame(int firstframe) {
  /* Loads in the first zones of data needed to do computation into the data buffer.
NOTE: Data buffer updated (by UpdateDataFrame()) before computations begin, 
  therefore each slide of data in buffer will get (virtually) shifted. Therefore, 
  this function does not completely fill buffer with all needed data. */
  
  int FirstFrameMod;
  char binfilename[LONGSTRING];
  
  FirstFrameMod = firstframe % Data_TRes;
  if(FirstFrameMod < 0) FirstFrameMod += Data_TRes;
  
  if(myrank == 0) cout << "Loading first Data Frame: " << firstframe << " with " << FirstFrameMod << flush;
  sprintf(binfilename,"%sf%04d.new",Data_Scratch,FirstFrameMod);
  
  ifstream Data_BinFileID(binfilename);
  if(!Data_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << binfilename << endl;
    exit(FAILURE);
  }
  
  Data_BinFileID.seekg(Data_HeaderSize,ios::beg);
  Data_BinFileID.read((char *)Data_Array[1], ND * sizeof(double)*Data_BlockSize);
  
  Data_BinFileID.close();
  if(myrank == 0) cout << " OK!" << endl;
}


void LoadFirstDataNestFrame(int firstframe) {
  /* Loads in the first zones of data needed to do computation into the data buffer.
NOTE: Data buffer updated (by UpdateDataFrame()) before computations begin, 
  therefore each slide of data in buffer will get (virtually) shifted. Therefore, 
  this function does not completely fill buffer with all needed data. */
  
  int FirstFrameMod;
  char binfilename[LONGSTRING];
  
  FirstFrameMod = firstframe % Data_TRes;
  if(FirstFrameMod < 0) FirstFrameMod += Data_TRes;
  
  if(myrank == 0) cout << "Loading first Data Frame: " << firstframe << " with " << FirstFrameMod << flush;
  
  for(int nn=0; nn<Nest_NumNests; nn++) {
   
    sprintf(binfilename,"%sf%04d_%04d.new",Data_Scratch,FirstFrameMod,Nest_List[nn]);
    
    ifstream Data_BinFileID(binfilename);
    if(!Data_BinFileID.is_open()) {
      cout << "\n Error opening file: \a" << binfilename << endl;
      exit(FAILURE);
    }
    
    int dblocksize;
    Data_BinFileID.read((char *)&DataNest_Min[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&DataNest_Max[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&DataNest_Res[1][nn*ND], ND * sizeof(int));
    Data_BinFileID.read((char *)&DataNest_Delta[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&dblocksize, sizeof(int));
    
 
    
    DataNest_Array[0][nn] = new double[ND*dblocksize];
    DataNest_Array[1][nn] = new double[ND*dblocksize];
    
    Data_BinFileID.read((char *)&DataNest_Array[1][nn][0], ND * dblocksize * sizeof(double));
    Data_BinFileID.close();
    
  }
  if(myrank == 0) cout << " OK!" << endl;
}


void UpdateDataNestFrame(int df2) {
  
  int df2mod = df2 % Data_TRes;
  if(df2mod < 0) df2mod += Data_TRes;
  
  char binfilename[LONGSTRING];
  
  if(myrank==0) cout << "Updating dataframe " << df2 << " with " << df2mod << flush;
  
  for( int nn=0; nn<Nest_NumNests ; ++nn) {
    
    sprintf(binfilename,"%sf%04d_%04d.new",Data_Scratch,df2mod,Nest_List[nn]);

    ifstream Data_BinFileID(binfilename);
    if(!Data_BinFileID.is_open()) {
      cout << "\n Error opening file: \a" << binfilename << endl;
      exit(FAILURE);
    }
    
    double *tempptr = DataNest_Array[0][nn];

    DataNest_Array[0][nn] = DataNest_Array[1][nn];
    DataNest_Array[1][nn] = tempptr;
    
    for(int d=0;d<ND;++d) {
      DataNest_Min[0][nn*ND+d] = DataNest_Min[1][nn*ND+d];
      DataNest_Max[0][nn*ND+d] = DataNest_Max[1][nn*ND+d];
      DataNest_Res[0][nn*ND+d] = DataNest_Res[1][nn*ND+d];
      DataNest_Delta[0][nn*ND+d] = DataNest_Delta[1][nn*ND+d];
    }

    int dblocksize;
    Data_BinFileID.read((char *)&DataNest_Min[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&DataNest_Max[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&DataNest_Res[1][nn*ND], ND * sizeof(int));
    Data_BinFileID.read((char *)&DataNest_Delta[1][nn*ND], ND * sizeof(double));
    Data_BinFileID.read((char *)&dblocksize, sizeof(int));
    Data_BinFileID.read((char *)&DataNest_Array[1][nn][0], ND * dblocksize * sizeof(double));
    Data_BinFileID.close();
    
  }
  
    
  if(myrank == 0) cout << " OK!" << endl;
  
}


void UpdateDataFrame(int df2) {
  
  int df2mod = df2 % Data_TRes;
  if(df2mod < 0) df2mod += Data_TRes;
  
  char binfilename[LONGSTRING];
  
  if(myrank==0) cout << "Updating dataframe " << df2 << " with " << df2mod << flush;
  
  sprintf(binfilename,"%sf%04d.new",Data_Scratch,df2mod);
  
  ifstream Data_BinFileID(binfilename);
  if(!Data_BinFileID.is_open()) {
    cout << "\n Error opening file: \a" << binfilename << endl;
    exit(FAILURE);
  }
  
  double *tempptr = Data_Array[0];
  
  Data_Array[0] = Data_Array[1];
  
  Data_Array[1] = tempptr;
  
  Data_BinFileID.seekg(Data_HeaderSize,ios::beg);
  Data_BinFileID.read((char *)Data_Array[1], sizeof(double)*Data_BlockSize * ND);
  
  Data_BinFileID.close();
  if(myrank == 0) cout << " OK!" << endl;
  
}


void AllocateMemoryForVelocityData() {
  /* Allocates memory for arrays that hold U, V, W data for 2 slides (i.e. the moving buffer), 
  and memory for a 3D array that serves as place holder during updates to the moving buffer. */ 
  
  if(myrank == 0) cout << "Allocating memory for velocity data ... ";
  
  if(Velocity_Format == 3) {
    DataNest_Array = new double**[2];
    if(DataNest_Array == NULL) {
      cout << "Memory Allocation Error" << endl;
      exit(1);
    }
    
    DataNest_Min = new double*[2];
    DataNest_Max = new double*[2];
    DataNest_Res = new int*[2];
    DataNest_Delta = new double*[2];
    
    for(int ss=0;ss<2;++ss) {
      DataNest_Array[ss] = new double*[Nest_NumNests];
      if(DataNest_Array[ss] == NULL) {
        cout << "Memory Allocation Error" << endl;
        exit(1);
      }
      
      DataNest_Min[ss] = new double[Nest_NumNests * ND];
      DataNest_Max[ss] = new double[Nest_NumNests * ND];
      DataNest_Res[ss] = new int[Nest_NumNests * ND];
      DataNest_Delta[ss] = new double[Nest_NumNests * ND];
    }
  }
  else {
    Data_Array = new double*[2];
    if(Data_Array == NULL) {
      cout << "Memory Allocation Error" << endl;
      exit(1);
    }
    for(int ss=0;ss<2;++ss) {
      Data_Array[ss] = new double[ND * Data_BlockSize];
      if(Data_Array[ss] == NULL) {
        cout << "Memory Allocation Error" << endl;
        exit(1);
      }
    }
  }
  
  
  
  if(myrank == 0) cout << "OK!" << endl;
}


void FreeMemoryForVelocityData() {
  
  if(myrank == 0) cout << "Freeing data memory ... " << flush; 
  if(Velocity_Format == 3) {
    
    for(int ss=0; ss<2; ++ss) {
      for(int nn=0; nn<Nest_NumNests; ++nn) 
        delete []DataNest_Array[ss][nn];
      delete []DataNest_Array[ss];
    }
    delete []DataNest_Array;
  }
  else {
    for(int ss=0; ss<2; ++ss) 
      delete []Data_Array[ss];
    
    delete []Data_Array;
  }
  if(myrank == 0) cout << "OK!" << endl;
  
}













