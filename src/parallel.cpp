// parallel.cpp
//
// parallel module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "parallel.h"

void MakeJobPackets( long int numpackets, long int numjobs, long int details[]) {
      
  long int numjobsperpkt = (long int)ceil((double)(numjobs)/(double)numpackets);
  long int npp = numpackets + numjobs - numpackets * numjobsperpkt;
  
  long int jobnumber = 0;
  
  for(int pktn=0; pktn<numpackets; ++pktn) {
    
    if(pktn < npp) 
      numjobs = numjobsperpkt;
    else 
      numjobs = numjobsperpkt - 1;
    
    details[pktn*3] = pktn;
    details[pktn*3+1] = jobnumber;
    details[pktn*3+2] = jobnumber+numjobs-1;
    jobnumber+=numjobs;
    
  }
}

