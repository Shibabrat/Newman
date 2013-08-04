// velocity.cpp
//
// velocity module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include <cmath>
#include <cstring>
using namespace std;
#include <gsl/gsl_matrix.h>
#include "velocity.h"
#include "boundary.h"
#include "data.h"
#include "parameters.h"

double *Plot_ArrayX;
double *Plot_ArrayU;
int     MyPlot_NumPoints;
ofstream  Plot_OutFileID;

void GetCodedVelocity(double t, double *X, double *dXdt) {
  
#if ND==2
   
    


  //Hill's Vortex
  /*dXdt[0] = 2.0*X[0]*X[1];
  dXdt[1] = -2.0*(2.0*X[0]*X[0]-1.0) - 2.0*X[1]*X[1];*/
  
  //Sigmoidal
  /*dXdt[0] = 0.4*tanh(10.0*(X[0]-0.5*t));
  dXdt[1] = 0.0;*/
  
  //Duffing Oscillator
  /*dXdt[0] = X[1];
  dXdt[1] = X[0] - X[0]*X[0]*X[0];*/
  
  //Campagnola
  /*dXdt[0] = 0.0;
   if(X[0]<Pendulum_Amplitude*t+5) 
    dXdt[1] = 1.0;
  else
    dXdt[1] = -1.0;*/
    
  //Pendulum
 
    
  dXdt[0] = X[1];
  dXdt[1] = -sin(X[0])-Pendulum_Amplitude*X[1]*sin(M_PI*t);
  
  //Hurricane

  /*const double sqrtterm = sqrt(1.0+4.0*Eye_alpha);
  const double beta = 1.0/sqrtterm;
  double yterm = X[1]-0.5*(1.0+sqrtterm);
  
  dXdt[0] = -yterm/(X[0]*X[0]+yterm*yterm+Eye_alpha) - beta;
  dXdt[1] = X[0]/(X[0]*X[0]+yterm*yterm+Eye_alpha) + Eye_Amplitude*X[1]*cos(Eye_Omega*t);*/

  //DNA
  
 /* const double epsilon = 1/1400.0;
  const double a = 7.0;
  const double x0 = 0.3; 
  const double factor = 2.0*epsilon*a/sqrt(DNA_N);

  vector<double> qhatvec(2*DNA_N,0.0);
  vector<double> qregvec(2*DNA_N,0.0);
  
  Lift(t,X,qhatvec);
  
  Hat2Reg(qhatvec,qregvec);
    
  dXdt[0]=X[1];
  dXdt[1]=0.0;
  for(int kk=0;kk<DNA_N;kk++) {
    double eterm = exp(-a*(1.0-cos(qregvec[kk])-x0));
    dXdt[1] += factor * (eterm-1.0) * eterm * sin(qregvec[kk]);
  }*/
  
  //Wedge
  
  /* complex<double> c1 = cos(W_lambda*X[1]);
  complex<double> s1 = sin(W_lambda*X[1]);
  
  complex<double> c2 = cos((W_lambda-2.0)*X[1]);
  complex<double> s2 = sin((W_lambda-2.0)*X[1]);
  
  double psi_theta = real(pow(X[0],W_lambda)*W_K0*(-W_alpha*W_lambda*s1 + (W_lambda-2.0)*W_beta*s2));
  double psi_r = real(W_lambda*pow(X[0],W_lambda-1.0)*W_K0*(W_alpha*c1 - W_beta*c2));
  
  dXdt[0] = 2.0*M_PI*psi_theta/X[0];
  dXdt[1] = -2.0*M_PI*psi_r/X[0];*/
  
  /*
   double weight = 0.1+exp(-(X[0]*X[0]+X[1]*X[1])/4.0);
   dXdt[0] = (-2*X[1]-X[1]*X[1])*weight;
   dXdt[1] = (2*X[0]-X[0]*X[1])*weight;*/
  
#elif ND==3
  
  //Hill Vortex
  /*dXdt[0] = 2.0*X[0]*X[2];
  dXdt[1] = 2.0*X[1]*X[2];
  dXdt[2] = -2.0*(2.0*(X[0]*X[0]+X[1]*X[1])-1.0)-2.0*X[2]*X[2];*/
  
 // Pendulum 3D 
 /* dXdt[0] = X[1];
  dXdt[1] = -sin(X[0])-0.65*X[1]*sin(M_PI*t);
  dXdt[2] = 0.0;*/
  
  /* ABC flow */
  const double Ap = sqrt(3);
  const double Bp = sqrt(2);

  dXdt[0] = (Ap+ABC_Amplitude*sin(M_PI*t))*sin(X[2]) + cos(X[1]);
  dXdt[1] = Bp*sin(X[0]) + (Ap+ABC_Amplitude*sin(M_PI*t))*cos(X[2]);
  dXdt[2] = sin(X[1]) + Bp*cos(X[0]);
  
  //Lorenz attractor
  /*dXdt[0] =10.0*(X[1]-X[0]);
  dXdt[1] = -X[0]*X[2]+28.0*X[0]-X[1];
  dXdt[2] = X[0]*X[1]-8.0/3.0*X[2]; */

#elif ND==4


	//Four body problem
    // given X[0], X[1], and t
    // FBP_mu, FBP_a, FBP_omega, FBP_theta0
	// possibly add FBP_mu3 parameter?
	double theta3 = FBP_omega * t + FBP_theta0;
	double X3 = FBP_a*cos(theta3);
	double Y3 = FBP_a*sin(theta3);
    
    double r1 = pow(((X[0] + FBP_mu)*(X[0] + FBP_mu) + (X[1]*X[1])), 1.5);
	double r2 = pow(((X[0] - (1-FBP_mu))*(X[0] - (1-FBP_mu)) + (X[1]*X[1])), 1.5);
	double r3 = pow(((X[0] - X3)*(X[0] - X3) + (X[1] - Y3)*(X[1] - Y3)),1.5) ;

    
    dXdt[0] = X[2];
    dXdt[1] = X[3];
	dXdt[2] = X[0] + 2.0*X[3] - ((1-FBP_mu)/r1)*(X[0] + FBP_mu) - (FBP_mu/r2)*(X[0] - (1-FBP_mu)) - (FBP_mu3/r3)*(X[0] - X3);
	dXdt[3] = X[1] - 2.0*X[2] - ((1-FBP_mu)/r1)*X[1] - (FBP_mu/r2)*X[1] - (FBP_mu3/r3)*(X[1] - Y3);


  // double den = 1.0/(1.0+TBP_ecc*cos(t));
  // double denA = pow(X[1]*X[1]+(X[0]+TBP_mu)*(X[0]+TBP_mu) , -1.5);
  // double denB = pow(X[1]*X[1]+(X[0]-1.0+TBP_mu)*(X[0]-1.0+TBP_mu) , -1.5);
  // dXdt[0] = X[2];
  // dXdt[1] = X[3];
  // dXdt[2] = ( X[0] - (1.0-TBP_mu)*(X[0]+TBP_mu)*denA - (X[0]-1.0+TBP_mu)*TBP_mu*denB ) * den + 2.0 * X[3];
  // dXdt[3] = ( X[1] - (1.0-TBP_mu)*X[1]*denA - X[1]*TBP_mu*denB ) * den - 2.0 * X[2];

#elif ND==6

  double den = 1.0/(1.0+TBP_ecc*cos(t));
  double denA = pow(X[2]*X[2]+X[1]*X[1]+(X[0]+TBP_mu)*(X[0]+TBP_mu) , -1.5);
  double denB = pow(X[2]*X[2]+X[1]*X[1]+(X[0]-1.0+TBP_mu)*(X[0]-1.0+TBP_mu) , -1.5);
  dXdt[0] = X[3];
  dXdt[1] = X[4];
  dXdt[2] = X[5];
  dXdt[3] = ( X[0] - (1.0-TBP_mu)*(X[0]+TBP_mu)*denA - (X[0]-1.0+TBP_mu)*TBP_mu*denB ) * den + 2.0 * X[4];
  dXdt[4] = ( X[1] - (1.0-TBP_mu)*X[1]*denA - X[1]*TBP_mu*denB ) * den - 2.0 * X[3];
  dXdt[5] = ( -X[2]*TBP_ecc*cos(t) - (1.0-TBP_mu)*X[2]*denA - X[2]*TBP_mu*denB ) * den;

#endif
  
}



void  Reg2Hat(vector<double>& qreg, vector<double>& qhat) {
  
  for(int jj=0;jj<DNA_N;jj++) {
    qhat[jj]=0.0;
    qhat[DNA_N+jj]=0.0;
    for(int ii=0;ii<DNA_N;ii++) {
      qhat[jj]+=(Utransform[jj][ii]*qreg[ii]);
      qhat[DNA_N+jj]+=(Utransform[jj][ii]*qreg[DNA_N+ii]);
    }
  }
}

void  Hat2Reg(vector<double>& qhat, vector<double>& qreg) {
  for(int jj=0;jj<DNA_N;jj++) {
    qreg[jj]=0.0;
    qreg[DNA_N+jj]=0.0;
    for(int ii=0;ii<DNA_N;ii++) {
      qreg[jj]+=(Utransform[ii][jj]*qhat[ii]);
      qreg[DNA_N+jj]+=(Utransform[ii][jj]*qhat[DNA_N+ii]);
    }
  }
}


void GenerateU() {
  
  /* Set DNA parameters */
    
  Utransform = new double*[DNA_N];
  for(int ii=0;ii<DNA_N;ii++) 
    Utransform[ii] = new double[DNA_N];
  
  double invsqrtN = 1.0/sqrt(DNA_N);
  double sqrt2N = sqrt(2.0)*invsqrtN;
  
  for(int k=0;k<DNA_N;k++)    /* w=0 */
    Utransform[0][k]=invsqrtN;
  
  for(int w=1;w<DNA_N/2;w++) {
    for(int k=0;k<DNA_N;k++) { 
      Utransform[w][k]=sqrt2N*cos(2.0*M_PI*k*w/DNA_N);
    }
  }
  
  for(int k=0;k<DNA_N;k++)     /* w=N/2 */
    Utransform[DNA_N/2][k]=invsqrtN*pow(-1.0,(double)k);
  
  for(int w=DNA_N/2+1;w<DNA_N;w++) {
    for(int k=0;k<DNA_N;k++) { 
      Utransform[w][k]=sqrt2N*sin(2.0*M_PI*k*w/DNA_N);
    }
  }
  
  DNA_q0.resize(2*DNA_N,0.0);
  
  /* Prescribe regular components */
    vector<double> qreg0(2*DNA_N,0.0);
    
    
    
    
    for(int nn=0;nn<DNA_N;++nn) {
      qreg0[nn] = 0.795398830184144/sqrt(DNA_N)+0.09*sin(2.0*M_PI*nn/DNA_N);
      qreg0[DNA_N+nn] = 0.0;
    }
        
    Reg2Hat(qreg0,DNA_q0);
    
    if(myrank==0) cout << "DNA_q0[0] = " << DNA_q0[0] << endl;
    
}


void GetNestedCartesianVelocity(double t,double *X,double *dXdt) {
  
  /* See which nest (if any) we are inside */
  
  int nest = 0;
  int inflag=0;
  
  for(int nn = Nest_NumNests; nn; --nn){
    inflag=0;
    nest = nn-1;
    for(int d=0;d<ND;++d) {
      if( (X[d] >= DataNest_Min[0][nest*ND + d]) && (X[d] <= DataNest_Max[0][nest*ND + d]) && (X[d] >= DataNest_Min[1][nest*ND + d]) && (X[d] <= DataNest_Max[1][nest*ND + d]) ) {
        inflag = 1;
      }
      else {
        inflag = 0;
        break;
      }  
    }
    if(inflag)
      break;
  }
  
  if(!inflag) {  
    /* We have no data for the requested point, Oh well, return 0. */ 
    for(int d=0;d<ND;++d) 
      dXdt[d]=0.0;
  }  
  else { 
    
    /* Use data from nest */
    
    //Define tloc  
    double tloc = (t-datatime1)/Data_TDelta*Time_Direction;
    double  uval[2][ND];
    
    for(int tf=0;tf<2;++tf) {     /* for each data frame */
      int ij0[ND];
      double loc[ND];
      GetNestedIJloc(X, loc, ij0, nest, tf);
     
      /*Linear interpolation by volumes */
      for(int v=0;v<Cell_NumVerts;++v) {
        Cell_Areas[v] = 1.0;
        for(int d=0;d<ND;++d)
          (v & Cell_Mask[d]) ? (Cell_Areas[v] *= loc[d]) : (Cell_Areas[v] *= (1.0-loc[d]));
      }
      
      int ij[ND];
      for(int dim=0;dim<ND;++dim) {  /* for u, v, w */
        
        uval[tf][dim] = 0.0;
        for(int v=0;v<Cell_NumVerts;++v) {
          for(int d=0;d<ND;++d) {
            ij[d] = ij0[d];
            if(v & Cell_Mask[d]) (ij[d])++;
          }
          
          uval[tf][dim] += Cell_Areas[v]*DataNest_Array[tf][nest][ij2f(ij,&DataNest_Res[tf][nest*ND],dim)];
          
        } // for each vertex
      } // for each velocity component 
    } // for each data frame  



for(int d=0;d<ND;++d) {
  if(Velocity_Null[d])
    dXdt[d] = 0.0;
  else
    dXdt[d] = uval[0][d]*(1.0-tloc) + uval[1][d]*tloc;
}




if(Atmos_Set) {
  dXdt[0]=dXdt[0]*180.0/(M_PI*Atmos_Radius*cos(X[1]*M_PI/180.0));
  dXdt[1]=dXdt[1]*180.0/(M_PI*Atmos_Radius);
}
  }  //else use nest data
}


void CopyPlottoWork(long int pktnum) {
  
  char syscommand[LONGSTRING];
  sprintf(syscommand,"cp %s%04ld.raw %s",Plot_Scratch,pktnum,Path_Work);
  system(syscommand);
  
}


void CleanUp_PlotVelocity() {
  
  Plot_OutFileID.close();
  
}

void Setup_PlotVelocity(long int *pktdetails) {
  
  GeneratePlotMesh(pktdetails);
  
  //char syscommand[LONGSTRING];
  //sprintf(syscommand,"rm -f %s*.raw",Plot_Scratch);
  
  char outfilename[LONGSTRING];
  sprintf(outfilename,"%s%04ld.raw",Plot_Scratch,pktdetails[0]);
  
  Plot_OutFileID.open(outfilename,ios::binary);
  if(!Plot_OutFileID.is_open()) {
    cout << "\n Error opening file: \a" <<  outfilename << endl;
    exit(FAILURE);
    
  }
}

/*#include <complex>
using namespace std;

const double W_K0 = 0.5/M_PI;
const double W_phi = 10.0*M_PI/180.0;
const complex<double> W_lambda(13.0795, 6.3844);
const complex<double> W_alpha = cos((W_lambda-2.0)*W_phi);
const complex<double> W_beta = cos(W_lambda*W_phi);*/


void  Lift(double t, double *q, vector<double>& qfull) {
  
    qfull[0]=q[0];
    qfull[DNA_N]=q[1];
    
    for(int ii=1;ii<DNA_N;ii++) {
      qfull[ii]=LinSol(t,ii);
      qfull[DNA_N+ii]=LinSol(t,DNA_N+ii);
    }
    
}

double  LinSol(double t,int w) {
  
  double  omega,sol;
  
  omega=sqrt(2.0*(1.0-cos(2.0*M_PI*(double)w/(double)DNA_N)));
  
  if(w==0) {
    sol = DNA_q0[0]+t*DNA_q0[DNA_N];
  }
  else if(w==DNA_N) {
    sol = DNA_q0[DNA_N];
  }
  else if(w<DNA_N) {
    sol = DNA_q0[w]*cos(omega*t)+DNA_q0[DNA_N+w]/omega*sin(omega*t);
  }
  else if(w>DNA_N) {
    sol = -DNA_q0[w-DNA_N]*omega*sin(omega*t)+DNA_q0[w]*cos(omega*t);
  }
  else {
    printf("Error in LinSol: w out of range\n");
    exit(1);
  }
  return(sol);
  
}




void GetCodedJacobian(double t, const double *X, double *dfdX, double *dfdt) {
  
#if ND==4
  double x = X[0];
  double y = X[1];
  double mu = TBP_mu;
  
  double fterm = 1.0/(1.0+TBP_ecc*cos(t));
  double denA = 1.0/(y*y+(1.0+x-mu)*(1.0+x-mu)); 
  double denA3 = denA*sqrt(denA);
  double denA5 = denA3*denA;
  double denB = 1.0/(y*y+(x-mu)*(x-mu)); 
  double denB3 = denB*sqrt(denB);
  double denB5 = denB3*denB;
  
  double A20 = ( 1.0 - mu*denA3 + 3.0*(1.0+x-mu)*(1.0+x-mu)*mu*denA5 + 3.0*(1.0-mu)*(x-mu)*(x-mu)*denB5 - (1-mu)*denB3 ) * fterm;
  double A31 = ( 1.0 - mu*denA3 + 3.0*y*y*mu*denA5 + 3.0*(1-mu)*y*y*denB5 - (1.0-mu)*denB3 ) * fterm;
  double A30 = ( 3.0*y*mu*(1.0+x-mu)*denA5 + 3.0*y*(1.0-mu)*(x-mu)*denB5 ) * fterm;
  double A21 = A30;
  
  gsl_matrix_view dfdX_mat = gsl_matrix_view_array (dfdX, ND, ND);
  gsl_matrix * m = &dfdX_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 1.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, 0.0);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 1.0);
  
  gsl_matrix_set (m, 2, 0, A20);
  gsl_matrix_set (m, 2, 1, A21);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, 2.0);
  
  gsl_matrix_set (m, 3, 0, A30);
  gsl_matrix_set (m, 3, 1, A31);
  gsl_matrix_set (m, 3, 2, -2.0);
  gsl_matrix_set (m, 3, 3, 0.0);
  
  double num = TBP_ecc*sin(t);
  
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = num * (x - (1.0+x-mu)*mu*denA3 - (1.0-mu)*(x-mu)*denB3  ) * fterm*fterm;
  dfdt[3] = num * (y - y*mu*denA3 - y*(1.0-mu)*denB3  ) * fterm*fterm;
#endif
  
  
}


void GetCartesianVelocity(double t,double *X,double *dXdt) {
	
  /* Simple test to see if we are outside the data domain */
  for(int d=0;d<ND;++d) {
    if(!Data_Periodic[d]){
      if( X[d] < Data_Min[d] || X[d] > Data_Max[d]) {
        for(int dim=0;dim<ND;++dim) {
          dXdt[dim] = 0.0;
        }
        return;
      }
    }
  }
  
  /* Test outside polygon */
  if(Boundary_Method == 2) {
    if(TestOutsidePolygon(X)) {
      for(int d=0;d<ND;++d) {
        dXdt[d] = 0.0;
      }
      return;
    }
  }
  
  //Define X index
  int ij0[ND];
  double loc[ND];
  
  GetIJloc(X, loc, ij0);
  
  /* Test outside staircase */
  if(Boundary_Method == 1) {
    if(Boundary_Mask[ij2f(ij0,Data_Res)]) {
      for(int d=0;d<ND;++d) {
        dXdt[d] = 0.0;
      }
      return;
    }
  }
  
  /*Linear interpolation by volumes */
  for(int v=0;v<Cell_NumVerts;++v) {
    Cell_Areas[v] = 1.0;
    for(int d=0;d<ND;++d)
      (v & Cell_Mask[d]) ? (Cell_Areas[v] *= loc[d]) : (Cell_Areas[v] *= (1.0-loc[d]));
  }
  
  //Define tloc  
  double tloc = (t-datatime1)/Data_TDelta*Time_Direction;
  double  uval[2][ND];
  int ij[ND];
  for(int dim=0;dim<ND;++dim) {
    for(int tf=0;tf<2;++tf) {
      uval[tf][dim] = 0.0;
      for(int v=0;v<Cell_NumVerts;++v) {
        for(int d=0;d<ND;++d) {
          ij[d] = ij0[d];
          if(v & Cell_Mask[d]) (ij[d])++;
        }
        uval[tf][dim] += Cell_Areas[v]*Data_Array[tf][ij2f(ij,Data_Res,dim)];
      }
    }
    
    if(Velocity_Null[dim])
      dXdt[dim] = 0.0;
    else
      dXdt[dim] = uval[0][dim]*(1.0-tloc) + uval[1][dim]*tloc;
    
  }
  
  if(Atmos_Set){
    dXdt[0]=dXdt[0]*180.0/(M_PI*Atmos_Radius*cos(X[1]*M_PI/180.0));
    dXdt[1]=dXdt[1]*180.0/(M_PI*Atmos_Radius);
  }
}  




void GeneratePlotMesh(long int *pktdetails) {
  
  if(myrank == 0) cout << "Generating Plot mesh ...";
  
  long int jbegin = pktdetails[1]; 
  long int jend = pktdetails[2];
  
  MyPlot_NumPoints = jend-jbegin+1;
  
  Plot_ArrayX = new double[ND*MyPlot_NumPoints];
  if(Plot_ArrayX == NULL) {
    cout << "Memory Allocation Error" << endl;
    exit(1);
  }
  
  int ij[ND];
  
  for(int pc=0;pc<MyPlot_NumPoints;++pc) {
    f2ij(jbegin+pc,ij,Plot_Res);
    for(int d=0;d<ND;++d) 
      Plot_ArrayX[ND*pc+d] = Plot_Min[d] + ij[d]*Plot_Delta[d];
    
  }
  
  Plot_ArrayU = new double[ND*MyPlot_NumPoints];
  
  if(myrank == 0) cout << "OK!" << endl;
}


void PlotVelocityFields(double t) {
  
  
  double X[ND];
  double U[ND];
  
  if(Plot_Velocity) {
    if(myrank == 0) 
      cout << "Plotting velocity field at t = " << t << endl; 
    
    for(int p = 0; p < MyPlot_NumPoints; ++p) {
      
      for(int d=0;d<ND;++d)
        X[d] = Plot_ArrayX[ND*p+d];
      
      GetVelocity(t, X, U);
      
      for(int d=0;d<ND;++d)
        Plot_ArrayU[d*MyPlot_NumPoints + p] = U[d];
      
    }
    
    for(int d=0;d<ND;++d){
      Plot_OutFileID.write((char *)&t,sizeof(double)); 
      Plot_OutFileID.write((char *)&MyPlot_NumPoints,sizeof(int)); 
      Plot_OutFileID.write((char *)&Plot_ArrayU[d*MyPlot_NumPoints],MyPlot_NumPoints*sizeof(double));  
    }
  }
  
  if(myrank == 0 || myrank == 1) {
    if(Query_Velocity) {
      GetVelocity(t, Query_X, U);
      int ij[ND];
      double loc[ND];
      if(Velocity_Format == 3) 
        GetNestedIJloc(Query_X, loc, ij, 0, 0);
      else
        GetIJloc(Query_X, loc, ij);
        
      cout << "\ntime = \t" << t << endl;
      double tloc = (t-datatime1)/Data_TDelta*Time_Direction;
      cout << "tloc = \t" << tloc << endl;
      
      cout << "ij = \t" << flush;
      for(int d=0;d<ND;++d)
        cout << ij[d] << " \t" <<flush;
      cout << endl;
      
      cout << "X = \t" << flush;
      for(int d=0;d<ND;++d) 
        cout << Query_X[d] << " \t" << flush;
      cout << endl;
     
      cout << "U = \t" << flush;
      for(int d=0;d<ND;++d) {
        cout << U[d] << " \t" << flush;
      }
      cout << endl;
    }
  }
}

void CombinePlotFiles(int numpackets) {
  
  char infilename[LONGSTRING];
  char systemcommand[LONGSTRING];
  
  double *U = new double[Plot_BlockSize];
  double frametime;
  double dataframetime;
  int     numplotpoints;
  
  cout << "Combining " << numpackets << " velocity plot files into raw binary format ... " << endl;
    
  ofstream fout(Plot_Output,ios::binary);
  if(!fout.is_open()) {
    cout << "\n Error opening file: \a" <<  Plot_Output << endl;
    exit(FAILURE);
  }
  
  int two = 2;
  fout.write((char *)&two, sizeof(int) );
  fout.write((char *)&Time_Origin[0],6*sizeof(int));
  fout.write((char *)&Output_TRes, sizeof(int) );
  fout.write((char *)&Atmos_Set, sizeof(int) );
  fout.write((char *)&Atmos_Radius, sizeof(double) );
  fout.write((char *)&Plot_BlockSize, sizeof(int) );
  fout.write((char *)&Plot_Res[0], ND*sizeof(int) );
  fout.write((char *)&Plot_Min[0], ND*sizeof(double) );
  fout.write((char *)&Plot_Max[0], ND*sizeof(double) );
  fout.write((char *)&Plot_Delta[0], ND*sizeof(double) );
  
  long int *fpos = new long int[numpackets];
  for(int p=0;p<numpackets;++p) fpos[p]=0;
  
  for(int tt=0;tt<Output_TRes;++tt) { 
    frametime = Output_T1 + tt*Output_TDelta;
    fout.write((char *)&frametime,sizeof(double));
    
    for(int d=0;d<ND;++d) {
      int count = 0;
      for(int p=0;p<numpackets;++p) {
        sprintf(infilename,"%s%04d.raw",Plot_Work,p);
        ifstream fin(infilename,ios::binary);
        if(!fin.is_open()) {
          cout << "\n Error opening file: \a" <<  infilename << endl;
          exit(FAILURE);
        }
        
        fin.seekg(fpos[p],ios::beg);
        
        fin.read((char *)&dataframetime, sizeof(double));
        fin.read((char *)&numplotpoints, sizeof(int));
        
        if(abs(dataframetime-frametime)>TINY_DATA) {
          cout << "\a WARNING: Frametimes in raw data do not agree: " << frametime << " " << dataframetime << endl;
        } 
        
        fin.read((char *)&U[count], numplotpoints * sizeof(double));
        
        count+=numplotpoints;
        fpos[p]=fin.tellg();
        fin.close(); 
        
      }
      if(count!=Plot_BlockSize) {
        cout << "Data sizes in raw data do not agree: " << count << " " << Plot_BlockSize << endl;
        exit(1);
      } 
      fout.write((char *)&U[0],Plot_BlockSize*sizeof(double));
    } 
  }
  
  fout.close();
  
  sprintf(systemcommand,"rm -f %s*.raw",Plot_Work);
  system(systemcommand);
  
  
  delete []fpos;
  cout << "Output file: " << Plot_Output << endl; 
}





