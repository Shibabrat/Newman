// boundary.cpp
//
// boundary module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "boundary.h"
#include "data.h"
#include <iostream>
#include <fstream>
using namespace std;

struct Boundary_List	*Boundary_Array;
int *Boundary_Mask;

void Setup_Boundary() {
    switch(Boundary_Method) {
        case 0 :
            break;
        case 1 : 
            GetBoundaryMask();
            break;
        case 2 : 
            ReadInBoundaryData(1);
            if(myrank==0) WriteOutTecBoundary();
            break;
            case 3 : //Analytical Boundary
            break;
            
            default : 
            cout << "ERROR: Incorrect Boundary_Method chosen:" << endl;
            cout << "0 no boundary, 1 staircase, 2 polygon, 3 analytical" << endl;
    }
}

void GetBoundaryMask() {
    
    if(myrank == 0) cout << "Creating Boundary Mask ... " << flush;
    char binfilename[LONGSTRING];
    int numcells = Data_BlockSize;
    
    double *veldata = new double[ND * Data_BlockSize];
    Boundary_Mask = new int[numcells];
    
    sprintf(binfilename,"%sf0000.new",Data_Work);
    
    ifstream Data_BinFileID(binfilename);
    if(!Data_BinFileID.is_open()) {
        cout << "\n Error opening file: \a" << binfilename << endl;
        exit(FAILURE);
    }
    
    Data_BinFileID.seekg(Data_HeaderSize,ios::beg);
    Data_BinFileID.read((char *)veldata, sizeof(double)*ND*Data_BlockSize);
    Data_BinFileID.close();
    
    int vertex;
    int vlim = 1;
    vlim <<= ND;
    
    int ij[ND];
    int ij0[ND];
    int mask[ND];
    for(int d=0;d<ND;++d) {
        mask[d] = 1;
        mask[d] <<= d;
    }
    /* For each cell */
    
    for(int nc = 0;nc<numcells;++nc) {
        
        int breakflag = 0;
        int edgeflag = 0; 
        Boundary_Mask[nc] = 0;
        
        /* Get the ij coordinates */
        f2ij(nc,ij0,Data_Res);
        
        for(int d=0;d<ND;++d) {
            if(ij0[d] > Data_Res[d]-2) {
                edgeflag = 1;
            }
        }
        if(edgeflag) {
            Boundary_Mask[nc] = 1;
            continue;
        }
        
        /* For each of the vertices on the cell */
        for(int v=0;v<vlim;++v) {
            
            for(int d=0;d<ND;++d) {
                ij[d] = ij0[d];
                if ( v & mask[d] ) (ij[d])++;
            }
            
            vertex = ij2f(ij,Data_Res);
            for(int d=0;d<ND;++d) {
                if(isnan(veldata[Data_BlockSize*d+vertex])) {
                    Boundary_Mask[nc] = 1;
                    breakflag = 1;
                    break;
                }
            }
            if(breakflag)
                break;
        }
    }
    
    delete []veldata;
    
    if(myrank == 0) cout << "OK!" << endl;
}

void CleanUpBoundary() {
    
    switch (Boundary_Method) {
            
        case 1 : 
            delete []Boundary_Mask;
            break;
            
        case 2 : 
            for(int list=0; list<Boundary_NumLists; ++list) {
                int numsegs = Boundary_Array[list].numsegments;
                for( int j=0; j<numsegs; ++j)
                    delete []Boundary_Array[list].X[j];
                delete []Boundary_Array[list].limits;
                delete []Boundary_Array[list].X;
            }
            delete []Boundary_Array;
            break;
    }
    
}

int TestOutsideAnalyticalBoundary(double *X, double t) {
    
#if ND==2  
    
    // DNA
    if(X[0]>M_PI || X[0]<0.0)
        return 1;
    
    // Wedge 
    /*double theta = atan2(X[1],X[0]);
     if(abs(theta)> 10.0*M_PI/180.0) return 1;*/
    
#elif ND==4
    
    /* Elliptic three-body problem */
    /*if (TBP_Ordinate == 0) {
     if(2.0*(X[2]-0.5*X[1]*X[1]+(0.5*X[0]*X[0]+0.5*X[3]*X[3]+0.5*TBP_mu-0.5*TBP_mu*TBP_mu+TBP_mu/sqrt(X[3]*X[3]+(X[0]-1.0+TBP_mu)*
     (X[0]-1.0+TBP_mu))+(1.0-TBP_mu)/sqrt(X[3]*X[3]+(X[0]+TBP_mu)*(X[0]+TBP_mu)))/(1.0+TBP_ecc*cos(t)))  < 0) 
     return 1;
     }
     
     if (TBP_Ordinate == 1) {
     if(2.0*(X[2]-0.5*X[1]*X[1]+(0.5*X[3]*X[3]+0.5*X[0]*X[0]+0.5*TBP_mu-0.5*TBP_mu*TBP_mu+TBP_mu/sqrt(X[0]*X[0]+(X[3]-1.0+TBP_mu)*(X[3]-1.0+TBP_mu))+(1.0-TBP_mu)/sqrt(X[0]*X[0]+(X[3]+TBP_mu)*(X[3]+TBP_mu)))/(1.0+TBP_ecc*cos(t))) < 0)
     return 1;
     }*/
    
    /* Four body problem */
    if(FBP_sunbdy2) {
        double distancefromsun2 = (X[0]+FBP_mu)*(X[0]+FBP_mu) + X[1]*X[1];
        if(distancefromsun2 < FBP_sunbdy2)
            return 1;
    }
    
    if(FBP_earthbdy2) {
        double distancefromearth2 = (X[0] - (1-FBP_mu))*(X[0] - (1-FBP_mu)) + (X[1]*X[1]);
        if(distancefromearth2 < FBP_earthbdy2)
            return 1;
    }
    
    if(FBP_moonbdy2) {
        double theta3 = FBP_omega * t + FBP_theta0;
        double X3 = FBP_a*cos(theta3);
        double Y3 = FBP_a*sin(theta3);
        double distancefrommoon2 = (X[0] - X3)*(X[0] - X3) + (X[1] - Y3)*(X[1] - Y3);
        if(distancefrommoon2 < FBP_moonbdy2)
            return 1;
    }
    
    
#elif ND==6
    if (TBP_Ordinate == 0) {
        double x = X[0];
        double vx = X[1];
        double E = X[2];
        double y = X[3];
        double z = X[4];
        double vz = X[5];
        double den = 1.0/(1.0+TBP_ecc*cos(t));
        double r1den = 1.0/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z);
        double r2den = 1.0/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z);
        if (2.0*(    E - 0.5*(vx*vx+vz*vz) + den*(0.5*(x*x+y*y)+0.5*TBP_mu*(1.0-TBP_mu)+TBP_mu*r2den+(1.0-TBP_mu)*r1den-0.5*TBP_ecc*z*z*cos(t)) ) < 0)
            return 1;
    }
    
    if (TBP_Ordinate == 1) {
        double y = X[0];
        double vy = X[1];
        double E = X[2];
        double x = X[3];
        double z = X[4];
        double vz = X[5];
        double den = 1.0/(1.0+TBP_ecc*cos(t));
        double r1den = 1.0/sqrt((x+TBP_mu)*(x+TBP_mu)+y*y+z*z);
        double r2den = 1.0/sqrt((x-1.0+TBP_mu)*(x-1.0+TBP_mu)+y*y+z*z);
        if (2.0*(    E - 0.5*(vy*vy+vz*vz) + den*(0.5*(x*x+y*y)+0.5*TBP_mu*(1.0-TBP_mu)+TBP_mu*r2den+(1.0-TBP_mu)*r1den-0.5*TBP_ecc*z*z*cos(t)) ) < 0)
            return 1;
    }
    
    
#endif
    
    return 0;
    
}


int TestOutsideStaircase(double *point) {
    int ij[ND];
    double xmod;
    for(int d=0;d<ND;++d) { 
        if(Data_Periodic[d]) {
            xmod = fmod( point[d]-Data_Min[d], Data_Period[d]);
            if(xmod<0) {
                xmod+= Data_Period[d];
            }
            xmod+=Data_Min[d];
        }
        else {
            xmod = point[d];
        }
        
        if(Data_NonUniformGrid[d]) {
            ij[d] = FindSlot(xmod,d,Data_Res);
        }  
        else {
            ij[d] = (int)floor((xmod-Data_Min[d])/Data_Delta[d]);
            
            if(ij[d] < 0) ij[d]=0;
            
            if(ij[d]>=(Data_Res[d]-1)){
                ij[d] =Data_Res[d]-2;
            }
        }
    }
    
    return( Boundary_Mask[ij2f(ij,Data_Res)] );
    
}

int TestOutsidePolygon(double *point) {
    /* Return 1 if outside, 0 if inside */
    
    if(QuickTestOutsidePolygon(0,point)) 
        return 1;
    
    if(LongTestOutsidePolygon(0,point)) 
        return 1;
    
    /* Check outside each polygon */
    if(Boundary_NumLists>1) 
        for(int j=1;j<Boundary_NumLists;++j) 
            if(!QuickTestOutsidePolygon(j,point)) 
                if(!LongTestOutsidePolygon(j,point)) 
                    return 1;
    
    return 0;
}

int QuickTestOutsidePolygon(int jbdy,double *point){
    
    for(int d=0;d<ND;++d) {
        if( point[d] < Boundary_Array[jbdy].limits[2*d] || point[d] > Boundary_Array[jbdy].limits[2*d+1]) {
            return(1);
        }
    }
    return 0;
}

int LongTestOutsidePolygon(int jbdy,double *point){
	
    /* Jul6 2007 For use with ND = 2 only */
    
    
	double	xnew,ynew;
	double	xold,yold;
	double	x1,y1;
	double	x2,y2;
	double	xt,yt;
	int		i;
	int		outside=1;
	int   numsegs;
    
    numsegs=Boundary_Array[jbdy].numsegments;
    
	xt=point[0];
	yt=point[1];
	
	xold=Boundary_Array[jbdy].X[numsegs-1][0];
	yold=Boundary_Array[jbdy].X[numsegs-1][1];
	
    
	for (i=0 ; i < numsegs ; i++) {
		xnew=Boundary_Array[jbdy].X[i][0];
		ynew=Boundary_Array[jbdy].X[i][1];
        
		if (xnew > xold) {
			x1=xold;
			x2=xnew;
			y1=yold;
			y2=ynew;
		}
		else {
			x1=xnew;
			x2=xold;
			y1=ynew;
			y2=yold;
		}
		if ((xnew < xt) == (xt <= xold)         /* edge "open" at left end */
			&& (yt-y1)*(x2-x1)
            < (y2-y1)*(xt-x1)) {
			outside=!outside;
		}
        xold=xnew;
        yold=ynew;
	}
    
    return(outside);
}

void  ReadInBoundaryData(int format) {
    
    switch (format) {
        case 1:
            ReadInTecBoundaryData();
            break;
        default:
            ReadInTecBoundaryData();
            break;
    }
    
}


void	ReadInTecBoundaryData() {
	
	char	buf[LONGSTRING];
    double min[ND], max[ND];
    for(int d=0; d<ND; ++d) {
        min[d]=HUGE_VAL;
        max[d]=-HUGE_VAL;
    }
    int   numsegs,listcount;
    
    
	if(myrank==0) cout << "Reading in Tecplot Boundary Data ... " << endl;
	
    Boundary_Array = new (nothrow) struct Boundary_List[MAXNUMISLANDS];
    if(Boundary_Array == NULL) 
        FatalError("Memory Allocation Error for Boundary: new failed.");
    
    
    /* Open boundary file */
    cout << Boundary_Input << endl;
    ifstream Boundary_InFileID(Boundary_Input);
    if( !Boundary_InFileID.is_open() ) {
		cout << "\nError opening boundary file: " << Boundary_Input << endl;
        exit(1);
	}
    
    /* Read in data */
    
    listcount=0;
    Boundary_InFileID.getline(buf,LONGSTRING,'\n');
    
    while(Boundary_InFileID.getline(buf,LONGSTRING,'\n') != NULL) {
        Boundary_InFileID.getline(buf,LONGSTRING,'\n');
        Boundary_InFileID.getline(buf,LONGSTRING,'\n');
        sscanf(buf,"%d\n",&numsegs);
        
        Boundary_Array[listcount].numsegments=numsegs;
        
        Boundary_Array[listcount].X = new (nothrow) double*[numsegs];
        if(Boundary_Array[listcount].X  == NULL) 
            FatalError("Memory Allocation Error for Boundary_Array[].X: new failed.");
        
        for(int jj=0;jj<numsegs;++jj) {
            Boundary_Array[listcount].X[jj] = new (nothrow) double[ND];
            if(Boundary_Array[listcount].X[jj] == NULL)
                FatalError("Memory Allocation Error for Boundary_Array[].X[]: new failed.");
        }
        
        Boundary_Array[listcount].limits = new (nothrow) double[2*ND];
        if(Boundary_Array[listcount].limits  == NULL) 
            FatalError("Memory Allocation Error for Boundary_Array.limits: new failed.");
        
        for(int kk=0;kk<numsegs;++kk) {
            
            for(int d=0;d<ND;++d)
                Boundary_InFileID >> Boundary_Array[listcount].X[kk][d];
            Boundary_InFileID.ignore(LONGSTRING,'\n');
            
            for(int d=0;d<ND;++d) {
                min[d]=fmin(min[d],Boundary_Array[listcount].X[kk][d]);
                max[d]=fmax(max[d],Boundary_Array[listcount].X[kk][d]);
            }
        }
        
        for(int d=0;d<ND;++d) {
            Boundary_Array[listcount].limits[2*d]=min[d];
            Boundary_Array[listcount].limits[2*d+1]=max[d];
        }
        ++listcount;
        if(myrank == 0) cout << "Read in " << numsegs << " segments" << endl;
    }
    
    Boundary_NumLists=listcount;
    if(myrank == 0) cout << "Number of BoundaryLists: " << listcount << endl;
	Boundary_InFileID.close();
    
}


void WriteOutTecBoundary() {
    /* Needs ND upgrade */ 
    
    char tecbdy_Output[100];
    
    if(strcmp(Path_Output,"pwd")) {
        sprintf(tecbdy_Output,"%sbdytec.dat",Path_Output);
    }
    else {
        sprintf(tecbdy_Output,"bdytec.dat");
    }
    
    
    ofstream tecID(tecbdy_Output);
    if(!tecID.is_open()) {
        cout << "Error opening boundary file: " << tecbdy_Output << endl;
        exit(1);
    }
    
    tecID << "VARIABLES=\"Longitude\"\"Latitude\"" << endl;
    for(int j=0;j<Boundary_NumLists;++j) {
        tecID << "GEOMETRY T=LINE, C=BLACK, CS=GRID, LT=0.2" << endl;
        tecID << "1" << endl;
        tecID << Boundary_Array[j].numsegments << endl;
        for(int k=0;k<Boundary_Array[j].numsegments;++k)
            tecID << Boundary_Array[j].X[k][0] << " " << Boundary_Array[j].X[k][1] << endl;
        
    }
    
    tecID.close();
    
}


