#include <cstdlib>
#include <cstdio>
#include "../source/ctqmcCHM.h"
#include "../source/ctqmcSIAM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  //double n = atof(argv[1]);
  double T = atof(argv[1]);

  // number of matsubara frequences to be used.
  // the fast fourier transform requires that this is a power of 2,
  // and that the number of tau points is the same
  int N = pow(2,13);

  //prepare ctqmcSIAM
  ctqmcSIAM siam;
  siam.M=50000000;
  siam.runParallel = true;
  siam.UseSmartNom = true; //this is by default
  siam.FreqCutoff = 12.0; //default 10.0

  char rfn[300];
  //sprintf(rfn,"/nfs/jaksa/QMC_triang/RUN%d",Nproc);
  //sprintf(rfn,"/nfs/jaksa/QMC_fixedDoping/run_ctqmc.n%.3f",n);
  sprintf(rfn,"/nfs/jaksa/QMC_fixedDoping/run_ctqmc.T%.3f",T);
  siam.runFolderName = rfn;
  siam.execName = "/nfs/jaksa/QMC_fixedDoping/ctqmc";

  // initialize the object that does the calculation
  ctqmcCHM chm;

  //assign siam to chm
  chm.ctqmcsiam = &siam;

  chm.SetBroydenOptions(false, false, 0); 
  chm.SetLoopOptions(35, 1.0e-4);

  chm.PatchTailWithAtomicLimit = false;
  chm.AtomicCutoff = 12.0;
  chm.PatchDelta = false;
  chm.SetUseBethe(true);
  chm.PHSymmetricCase = false;

  chm.UseFixedMuSIAMRun = false; //if true set mu, if false set n and initial guess mu
  chm.PrintMuHistory = true;
  chm.UseBroydenForMu = true;
  chm.BroydenStartDiff = 5e-2;
  chm.UseBroydenForMu = true;
  chm.UseSmart_c = true;
  chm.UseSmart_M = true;
  chm.Smart_it = 20;
  chm.Smart_smallM = 5000000;
  chm.Smart_largeM = 50000000; 


  chm.UseLambdaCalculator = false;
  chm.LC->SetDoOutput(false);

  chm.PrintIntermediate = 0;
  //chm.IntermediateCutoff = 10.0; //default 10.0
  chm.PrintDiffs = true; 

  chm.UseIPT = false;

  // iterate over parameters
  //for(double T=0.01; T<0.65; T*=2.0)
  for(double n=0.5; n<0.7; n+=0.02)
  {
    // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
    IAGRID iagrid(N,N, T);
    // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
    IAResult iaresult(&iagrid);
    // initialize FFT. this is also dependent on temperature
    //chm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); //not needed for ctqmc

    char muFN[300];
    //sprintf(muFN,"mu_vs_n.n%.3f",n);
    sprintf(muFN,"mu_vs_n.T%.3f",T);         
    FILE* muFile = fopen(muFN,"w");
    fclose(muFile);
 
    //for(double n=0.5; n<0.7; n+=0.02)
    iaresult.mu = 0.75 + 3.0*(n-0.5); //initial guess for mu
    for(double U=1.5; U<4.01; U+=0.5)
    {       
      char bareFN[300];
      sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
      char FN[300];
      sprintf(FN, "ctqmcCHM%s",bareFN);     
      chm.intermediateName = FN;    

      chm.SetParams(U, T, 0.5);
      iaresult.n = n;

      //--------- RUN THE CODE---------//
      printf("about to run CTQMC\n");
      chm.Run(&iaresult);

      iaresult.PrintMinimal(FN,1.0e+300);

      muFile = fopen(muFN,"a");
      fprintf(muFile,"%.15le %.15le %.15le %.15le\n", n, U, iaresult.mu, iaresult.n0); //here n0 is the actual n, and in IPT it is n(G_0)
      fclose(muFile);

    }
  }
}

