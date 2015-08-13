#include "ctqmcDHM.h"
#include "ctqmcSIAM.h"
#include "IAResult.h"
#include "IAresArray.h"
#include "routines.h"

#include <cstdlib>
#include <cstdio>

ctqmcDHM::ctqmcDHM(): IADHM()
{
  UseSeparateFolders = true;
  UseIPT = false;
}


ctqmcDHM::~ctqmcDHM()
{ //was eva
} 

bool ctqmcDHM::SolveSIAM()
{
  if (UseIPT) return IADHM::SolveSIAM();

  if (ctqmcsiam==NULL) 
  { printf("ctqmcsiam not initialized! exiting to system\n");
    exit(0);
  }
  string original_rfn = ctqmcsiam->runFolderName;

  bool err = false;
  for(int id=0; id<Nsites; id++)
  {  ctqmcsiam->U = U;
     ctqmcsiam->beta = 1.0/T;
     ctqmcsiam->mu = a->r[id].mu;
     ctqmcsiam->PatchTailWithAtomicLimit = PatchTailWithAtomicLimit;
     ctqmcsiam->AtomicCutoff = AtomicCutoff;

     if (UseSeparateFolders)
     { char rfn[300];
       sprintf(rfn,"%s/run_ctqmc.%d",ctqmcsiam->runFolderName.c_str(), id);
       ctqmcsiam->runFolderName = rfn;
     }

     err = ctqmcsiam->Run(&(a->r[id]));
 
     a->r[id].n = ctqmcsiam->GetOccupancy();

     ctqmcsiam->runFolderName = original_rfn;
  }

  return err;

}
