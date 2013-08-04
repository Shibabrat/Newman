#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include  <mpi.h>
#include <map>
#include <ctime>

#include "macros.h"
#include "boundary.h"
#include "tracers.h"
#include "integrate.h"
#include "velocity.h"
#include "errors.h"
#include "ftle.h"
#include "parameters.h"
#include "data.h"
#include "globals.h"
#include  "parallel.h"
#include "mapcoord.h"

void  DoPacket(long int *pktdetails);

int  main(int argc, char *argv[]) {
  
  /*Initialize MPI */
  MPI_Status status;	
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  /* Usage instructions */
  if((myrank == 0) & (argc != 2)) {
    cout << "newmanND usage: >> newmanND parametersfile.in" << endl;
    exit(1);
  }
  
  /* Set Parameters */
  SetDefaultParameters();
  ReadInParameters(argv);
  SetDerivedParameters();
  
  /* Convert ASCII Data to binary */
  if(myrank == 0) 
    CheckParameters(); 
  
  if(Velocity_Format > 1) {
    if(Data_Format) {
      ReadInASCIIVelocityData(Data_Format,NumProcs); 
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  long int NumPackets = 1;  
  /********** MASTER ***************************/ 
  if(myrank==0) {
    
    /* Vital Statistics */
    long int NumJobs = 1;
    if(FTLE_Compute) NumJobs = FTLE_BlockSize*Output_TRes;
    if(Trace_Compute) NumJobs = Trace_NumDrifters;
    if(Plot_Velocity) NumJobs = Plot_BlockSize;
    
    NumPackets = (long int) ((NumProcs-1)*Parallel_LoadRatio);
    if(NumPackets==0) 
      NumPackets = 1;
    cout << "Running with MPI on " << NumProcs << " processes." << endl;
    cout << "Dimensions = " << ND << endl; 
    cout << "NumJobs = " << NumJobs << endl;
    cout << "NumPackets = " << NumPackets << endl;
    
    /* Make up job packets */
    long int *PacketInfo = new long int[3*NumPackets];   // Each packet has [ PacketID, jbegin, jend ]
   
    MakeJobPackets(NumPackets,NumJobs,PacketInfo);
    
    long int KillSignal=-1; 
    
    /*******************************************/
    /* SINGLE PROCESSOR */
    /*******************************************/
    if(NumProcs == 1) {  /*Ooops.  No slaves, Master better jump in.*/
      
      DoPacket(PacketInfo);
      
    }
    else {
      
      /***********************************/
      /* MULTIPLE PROCESSORS */
      /***********************************/
      
      /* Initial assignment of jobs to each processor */
      int PacketsSent=0;
      while(PacketsSent < (fmin(NumPackets,NumProcs-1)) ) {
        MPI_Send(&PacketInfo[3*PacketsSent],3,MPI_LONG,PacketsSent+1,1,MPI_COMM_WORLD);      
        cout << "Dispatched packet " << PacketsSent+1 << " to process " << PacketsSent+1 << endl; 
        ++PacketsSent;
      }
      
      // Listen, Record results, and send out more jobs
      
      for(int PacketsRecvd = 0;PacketsRecvd<NumPackets;) {
        
        int PacketStatus;
        MPI_Recv(&PacketStatus,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        int PacketID = status.MPI_TAG;
        int ProcID = status.MPI_SOURCE;
        
        
        ++PacketsRecvd;
        cout << PacketsRecvd << " Packet " <<  PacketID << " on process " << ProcID << " is finished: \t" << PacketStatus << " s" << endl;
        
        if(PacketsSent < NumPackets) {
          MPI_Send(&PacketInfo[3*PacketsSent],3,MPI_LONG,ProcID,1,MPI_COMM_WORLD);
          cout << "Dispatched packet " << PacketsSent+1 << " to process " << ProcID << endl; 
          ++PacketsSent;
        }
        else {
          MPI_Send(&KillSignal,1,MPI_LONG,ProcID,1,MPI_COMM_WORLD);
        }
      }
    }
  }
  else {   
    /*********************SLAVE NODES***********************************/
    
    /* Start loop to get jobs from MASTER */
    while(1) {
      
      /* get details for the job to complete */
      long int PacketDetails[3];
      MPI_Recv(&PacketDetails,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      long int PacketID = PacketDetails[0];
 
      /* Kill signal? */
      if(PacketID<0) break;
      
      /* Start running the packet */
      int t_on = clock();
      DoPacket(PacketDetails);
      int t_off = clock();
      
      /* Tell Master Node that the packet is completely finished */
      int PacketStatus = (t_off - t_on)/CLOCKS_PER_SEC;
      MPI_Send(&PacketStatus,1,MPI_INT,0,PacketID,MPI_COMM_WORLD);
      
    } /* End While */
} /* End Slave Node */

MPI_Barrier(MPI_COMM_WORLD);
    if(FTLE_Compute) {
      
      
      for(int ss=myrank;ss<Output_TRes;ss+=NumProcs) {
        CombineFTLEFiles(ss);
        OutputFTLEtoFile(ss);
      }
      
    }

    if(Trace_Compute) 
      if(myrank==0) 
        CombineTraceFiles(NumPackets);

    if(Plot_Velocity) 
      if(myrank==0)
        CombinePlotFiles(NumPackets); 

if(myrank == 0) cout << "Finalizing ... " << endl;

MPI_Finalize();

if(myrank == 0) cout << "DONE. \a" << endl;

return 0;

}




void DoPacket(long int *pktdetails) {

  double t, t1, t2;
  
  double TempPoint[ND];                   /* Used to save integrated position of FTLE point */

  map<int,struct FTLE_Drift> TempFTLE_Dfts;
  
  long int packetID = pktdetails[0];
  long int jbegin = pktdetails[1]; 
  long int jend = pktdetails[2];
  
  //GenerateU();
  
  SetMyDerivedParameters(pktdetails);
  
  ReadInNonUniformGrids();
  
  if(Velocity_Format > 1) {
    (Velocity_Format == 3) ? ReadInDataNestParameters() : ReadInDataParameters();
    AllocateMemoryForVelocityData();
  }
  
  if(Boundary_Method) Setup_Boundary();
  
  if(Plot_Velocity) Setup_PlotVelocity(pktdetails);
  
  if(FTLE_Compute) 
    AllocateMemoryForFTLE(pktdetails);
  
  if(FTLE_Compute || Trace_Compute) SetUpIntegrator();
  
  if(Velocity_Format > 1) {
    
    Velocity_Format == 3 ? LoadFirstDataNestFrame(Data_FirstFrame) : LoadFirstDataFrame(Data_FirstFrame);
    
    if(Trace_Compute) Setup_Trace(pktdetails);
    
    /* For each frame of data */
    if(myrank == 0) {
      cout << "DFF " << Data_FirstFrame << endl;
      cout << "DLF " << Data_LastFrame << endl;
    }
    for(int df1 = Data_FirstFrame; (Data_LastFrame-df1)*Time_Direction >= 0; df1+=Time_Direction) { 
      if(myrank==0) {
        double completepct = 100.0 * (double)(df1 - Data_FirstFrame) / (double)(Data_LastFrame - Data_FirstFrame + Time_Direction);
        cout << setprecision(4) << fixed << completepct << "% Complete" << endl;
      }
      Velocity_Format == 3 ? UpdateDataNestFrame(df1+Time_Direction) : UpdateDataFrame(df1+Time_Direction);
      
      datatime1 = Data_TMin + df1 * Data_TDelta;
      datatime2 = datatime1 + Time_Direction * Data_TDelta;
      
      if(myrank == 0) cout << setprecision(4) << showpoint << datatime1 << " to " << datatime2 << endl;
      
      /************* FTLE_COMPUTE ****************/
      if(FTLE_Compute) {
        long int ssbegin = (long int)(jbegin/FTLE_BlockSize);
        long int ssend = (long int)(jend/FTLE_BlockSize);
        
        /* For each output time slide */
        for(int ss = 0; ss < MyOutput_TRes; ++ss) { 
          
          /* If slide not complete */ 
          if(FTLE[ss].completestatus) {
            
            long int j1, j2;
            ss==0 ? (j1 = jbegin-(ssbegin*FTLE_BlockSize)) : (j1=0);
            ss==(MyOutput_TRes-1) ? (j2 = jend-(ssend*FTLE_BlockSize)) : (j2=FTLE_BlockSize-1);
            
            /* Consider only slides launched, or needing to be launched */ 
            if( (datatime2 - FTLE[ss].launchtime)*Time_Direction > 0.0) {
              
              /* Get slide data (FTLE_Array) from binary file */
              ReadInFTLESlide(ss,packetID);
              
              /* Initialize slide if it hasn't been launched yet */
              if(FTLE[ss].completestatus == 2) {
                FTLE[ss].completestatus = 1;
                t1 = FTLE[ss].launchtime;
                if(myrank == 0) cout << "Initiated slide " << ss << " for launch at time " << t1 << endl;
              }
              else {
                t1 = datatime1;
              }
              
              /* Integrate Slide */
              if(myrank==0) cout << "Integrating slide " << ss << " " << flush;
              
              t2 = Time_Direction * fmin(Time_Direction*datatime2, Time_Direction*FTLE[ss].stoptime);
              
              int statusstep = FTLE_Dfts.size()/10;
              int count = 0;
              
              map<int,struct FTLE_Drift>::iterator it;
              for(it=FTLE_Dfts.begin() ; it != FTLE_Dfts.end(); it++) {
                
                ++count;
                int f = (*it).first;
                
                if(myrank == 0) if(!(count % statusstep)) cout << "." << flush;
                
                /* Push the drifter forward */
                for(int d=0;d<ND;++d)
                  TempPoint[d] = FTLE_Dfts[f].X[d];
                
                Integrate(TempPoint, t1, t2);
                
                /* Drifter left domain */
                if(TestOutsideDomain(TempPoint,t2)) { 
                  
                  GetFTLE(f, FTLE_IntTLength, j1, j2);  /* get FTLE for all FTLE_Pts that have f as a neighbor*/
                  
                } 
                else {
                  /* Add drifter to new list of drifters */
                  for(int d=0;d<ND;++d) 
                    TempFTLE_Dfts[f].X[d] = TempPoint[d];
                }
              } // for each drifter
              if(myrank == 0) cout << endl;
              
              /* Copy contents of FTLE_NewArray to FTLE_Array */ 
              
              FTLE_Dfts=TempFTLE_Dfts;
              
              /* Write contents of FTLE_Array to binary file for slide ss */
              
              WriteOutFTLESlide(ss,packetID); 
            } //if integrate slide
            
            if((FTLE[ss].stoptime - datatime2)*Time_Direction < TINY_FTLE) {  
              
              FTLE[ss].completestatus = 0;
              if(myrank == 0) cout << "Terminating slide " << ss << endl;
              ComputeFTLE(ss,j1,packetID);
            } /* If terminate slide */
          } /* If slide not complete */
        } /* For each output time slide */
        
        
      } /* If FTLE_Compute */
      
      /************* TRACE_COMPUTE ***************/
      if(Trace_Compute) {
        t1 = datatime1;
            
        /* Output if necessary */
        if(t1 == nextoutputtime) {
          OutputTracers(t1);
          nextoutputtime += (Output_TDelta*Time_Direction);
        }
        
        /*Integrate up to next event*/
        t2 = Time_Direction*fmin(Time_Direction*datatime2, Time_Direction*Int_T2);
        t  = Time_Direction*fmin(Time_Direction*t2,Time_Direction*nextoutputtime);
        
        while((t2-t1)*Time_Direction > 0) {
          
          IntegrateTracers(t1, t);
          t1 = t;
          
          if( fabs(nextoutputtime-t1) < TINY_DATA) {
            OutputTracers(t1);
            nextoutputtime += (Output_TDelta*Time_Direction);
          }
          t  = Time_Direction*fmin(Time_Direction*t2,Time_Direction*nextoutputtime);
        }
        
      }
      
      /************* PLOT_VELOCITY ****************/ 
      if(Plot_Velocity) {
        
        t1 = datatime1;
        
        /* Output if necessary */
        if(t1 == nextoutputtime) {
          PlotVelocityFields(t1);
          nextoutputtime += (Output_TDelta*Time_Direction);
        }
        
        /*Find up to next event*/
        t2 = Time_Direction*fmin(Time_Direction*datatime2, Time_Direction*Int_T2);
        t  = Time_Direction*fmin(Time_Direction*t2,Time_Direction*nextoutputtime);
        
        while((t2-t1)*Time_Direction > 0) {
          t1 = t;
          if( fabs(nextoutputtime-t1) < TINY_DATA) {
            PlotVelocityFields(t1);
            nextoutputtime += (Output_TDelta*Time_Direction);            
          }
          t  = Time_Direction*fmin(Time_Direction*t2,Time_Direction*nextoutputtime);
        }
      }
    }
    
    /***** CLEAN_UP *******/
    
    FreeMemoryForVelocityData(); 
    
  }
  else {  // Use Analytical velocity
    
  
    /************* FTLE_COMPUTE ****************/
    if(FTLE_Compute) {
      
      /* For each output time slide */
      for(int ss = 0; ss < MyOutput_TRes; ss++) { 
        
        long int ssbegin = (long int)(jbegin/FTLE_BlockSize);
        long int j1;
        ss==0 ? (j1 = jbegin-(ssbegin*FTLE_BlockSize)) : (j1=0);
        
        /* Get slide data (FTLE_Array) from binary file */
        ReadInFTLESlide(ss,packetID);
        
        /* Initialize slide if it hasn't been launched yet */
        
        t1 = FTLE[ss].launchtime;
        t2 = FTLE[ss].stoptime;
        if(myrank == 0) cout << "Initiated slide " << ss << " for launch at time " << t1 << " " << flush;
        
        int statusstep = FTLE_BlockSize/10;
        int count=0;
        map<int,struct FTLE_Drift>::iterator it;
        for(it=FTLE_Dfts.begin() ; it != FTLE_Dfts.end(); it++) {
          ++count;
          int f = (*it).first;
          
          if(myrank == 0) {
            if(!(count % statusstep)) cout << "." << flush;
          }
          
          Integrate(FTLE_Dfts[f].X , t1, t2);
        } 
        ComputeFTLE(ss,j1,packetID);
        if(myrank == 0) cout << endl;
        
      }
    }
    
    
    /************* TRACE_COMPUTE ***************/
    if(Trace_Compute) {
      if(Trace_Compute) Setup_Trace(pktdetails);
      OutputTracers(Int_T1);
      t1 = Int_T1;
      t2 = Int_T1 + Output_TDelta*Time_Direction;
      for(int tt=0;tt<Output_TRes-1;++tt) {
        if(myrank==0) cout << "Integrating trajectories from " << t1 << " to " << t2 << endl;
        
        IntegrateTracers(t1, t2);
        
        OutputTracers(t2);
        t1=t2;
        t2+=(Output_TDelta*Time_Direction);
      }
    }
    
    /************* PLOT_VELOCITY ****************/ 
    if(Plot_Velocity) {
      nextoutputtime = Int_T1;
      for(int tt=0; tt<Output_TRes; ++tt) {
        PlotVelocityFields(nextoutputtime);
        nextoutputtime += (Output_TDelta*Time_Direction);
      }
    }
  }
  
  if(Plot_Velocity) {
    CleanUp_PlotVelocity();
    CopyPlottoWork(packetID);
  }
  if(Trace_Compute) {
    CleanUpTrace();
    CopyTracetoWork(packetID);
  }
  
  if(Boundary_Method) CleanUpBoundary();
  
  if(FTLE_Compute || Trace_Compute) CleanUpIntegrator();
  
  if(FTLE_Compute) {
    CopyFTLEtoWork(packetID);
    FreeMemoryForFTLE();  
  }
}




