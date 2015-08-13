#include "ctqmcCHM.h"
#include "ctqmcSIAM.h"
#include "IAResult.h"
#include "routines.h"
//#include "Broyden.h"
#include <cstdlib>
#include <cstdio>

ctqmcCHM::ctqmcCHM(): IACHM()
{
  c = 1.0; //default value

  PrintMuHistory = true;
  UseBroydenForMu = true;
  BroydenStartDiff = 3e-2;
  BroydenStatus = 0;

  UseIPT = false;

  LastSign=0;
  UseSmart_c = false;
  UseSmart_M = false;
  Smart_it = 10;
  Smart_smallM = 5e+6;
  Smart_largeM = 50e+6; 
}


ctqmcCHM::~ctqmcCHM()
{ //was eva
} 



bool ctqmcCHM::SolveSIAM()
{
  if (UseIPT) return IACHM::SolveSIAM();

  char muHistoryFN[300];
  sprintf(muHistoryFN,"mu_history.n%.3f.U%.3f.T%.3f",r->n,U,T);

  if (Iteration==1) LastSign=0;

  if ((Iteration==1)and(PrintMuHistory))
  {
    FILE* f= fopen(muHistoryFN,"w");
    fclose(f);
  }

  if (ctqmcsiam==NULL) 
  { printf("ctqmcsiam not initialized! exiting to system\n");
    exit(0);
  }
  if (UseSmart_M)
  { if(Iteration<Smart_it)
      ctqmcsiam->M = Smart_smallM;
    else
      ctqmcsiam->M = Smart_largeM;
  }
  ctqmcsiam->U = U;
  ctqmcsiam->beta = 1.0/T;
  ctqmcsiam->mu = r->mu;
  ctqmcsiam->PatchTailWithAtomicLimit = PatchTailWithAtomicLimit;
  ctqmcsiam->AtomicCutoff = AtomicCutoff;

  bool err = ctqmcsiam->Run(r);
 
  if (UseFixedMuSIAMRun)
    r->n = ctqmcsiam->GetOccupancy();
  else
  { double n = ctqmcsiam->GetOccupancy();
    
    r->n0 = n; // save the actual n from calculation in n0, altough it can be found in observables as well. not to be confused with n(G_0) in IPT

    if ( (abs(r->n - n) < BroydenStartDiff) and (BroydenStatus==0) and ((Iteration>Smart_it) or (not UseSmart_M)) )
    {  B.SetParameters(1, MAX_ITS, 1.0, 0.01, 5e-5);
       B.TurnOn(Iteration-1);
       BroydenStatus++;
    }
    if ((UseSmart_c)and(BroydenStatus==0))
    { if (LastSign==0) 
      { c=0.2;
        LastSign = int_sign(r->n - n);
      }
      else
        if (LastSign==int_sign(r->n - n)) c*=1.5;
        else LastSign *= -1;
    }
    
    r->mu += c*(r->n - n); 
    
    if ( (UseBroydenForMu)and(BroydenStatus==1) )
    { complex<double> V[1];
      V[0] = r->mu;

      int conv = B.CalculateNew(V, Iteration);
      if (conv==1)
      { printf("mu converged, turning off broyden");
        BroydenStatus++;
      }
      if (abs(real(V[0]))<7.0)
        r->mu = real(V[0]);
      else 
      { printf("Broyden wandered off, turning off broyden");
        BroydenStatus++;
      }
    }

    if (PrintMuHistory)
    { FILE* f= fopen(muHistoryFN,"a");
      fprintf(f,"%d %.15le %.15le\n", Iteration, r->mu, n);
      fclose(f);
    }
  }
  
  return err;

}
