#include "../source/ctqmcCHM.h"
#include "../source/ctqmcSIAM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char* argv[])
{
  //if(argc<2) exit(0);
  //int Nproc = atof(argv[1]);


  // number of matsubara frequences to be used.
  // the fast fourier transform requires that this is a power of 2,
  // and that the number of tau points is the same
  int N = pow(2,13);


  //prepare ctqmcSIAM
  ctqmcSIAM siam;
  siam.M=20000000;
  siam.runParallel = true;

  char rfn[300];
  //sprintf(rfn,"/nfs/jaksa/QMC_triang/RUN%d",Nproc);
  sprintf(rfn,"/nfs/jaksa/QMC_triang/run_ctqmc");
  siam.runFolderName = rfn;
  siam.execName = "/nfs/jaksa/QMC_triang/ctqmc";


  // initialize the object that does the calculation
  ctqmcCHM chm;

  //assign siam to chm
  chm.ctqmcsiam = &siam;

  // set DMFT loop option
  chm.SetBroydenOptions(false, false, 0); 
  //iachm.SetBroydenOptions(true, false, 1.0e-4); 
  chm.SetLoopOptions(30, 1.0e-4);
  chm.PatchTailWithAtomicLimit = false;
  chm.AtomicCutoff = 12.0;
  chm.SetUseBethe(true);
  chm.UseFixedMuSIAMRun = false; //if true set mu, if false set n and initial guess mu


  chm.SetPrintOutOptions(false, false); // sth, halt on iterations

  // iterate over parameters
  for(double T=0.05001; T<0.1; T+=1000.05)
  {
    /*siam.nom = (int) ( 0.5 * ( 10.0/(pi*T) - 1.0 ) );
    if (siam.nom == 0) siam.nom = 1;
    if (siam.nom > 100) siam.nom = 100;*/ //this is now handled by UseSmartNom in ctqmcSIAM. Set also FreqCutoff



    // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
    IAGRID iagrid(N,N, T);
    // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
    IAResult iaresult(&iagrid);
    // initialize FFT. this is also dependent on temperature
    chm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); 


    double n= 0.6;
    for(double U=2.0; U<4.0; U+=10000.5)
    { char bareFN[300];
      sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
 
      chm.SetParams(U, T, 0.5);
    
      //--------- RUN THE CODE---------//
      iaresult.mu = U/2.0;
      iaresult.n = n;
      chm.Run(&iaresult);

      // printout result
      char FN[300];
      sprintf(FN, "ctqmcCHM%s",bareFN);
      iaresult.PrintMinimal(FN);
        
        // use Pade to analytically continue Green's funtion to the real axis
        //sprintf(FN, "Gw%s.IPT",bareFN);
        //PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );
    }
  }


}


