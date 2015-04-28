#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include <cstdio>

int main()
{

  // number of matsubara frequences to be used.
  // the fast fourier transform requires that this is a power of 2,
  // and that the number of tau points is the same
  int N = pow(2,13);

  // initialize the object that does the calculation
  IACHM iachm;

  // set DMFT loop option
  iachm.SetBroydenOptions(false, false, 0); 
  //iachm.SetBroydenOptions(true, false, 1.0e-4); 
  iachm.SetLoopOptions(200, 1.0e-9);
  iachm.PatchTailWithAtomicLimit = false;
  iachm.AtomicCutoff = 12.0;

  iachm.SetPrintOutOptions(false, false); // sth, halt on iterations


  // iterate over parameters
  for(double T=0.2001; T>0.0; T-=0.05)
  {
    // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
    IAGRID iagrid(N,N, T);
    // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
    IAResult iaresult(&iagrid);
    // initialize FFT. this is also dependent on temperature
    iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); 
    

    for(double n=0.5; n<0.7; n+=1000.05)
    {
      // if non interacting DOS is symetric around w=0 and n=0.5 then the calculation is simplified. use particle hole symmetry
      if (n!=0.5) iachm.PHSymmetricCase = false;
      else 
      { iachm.PHSymmetricCase = true;
        iachm.SetForceSymmetry(false);
      }
      iaresult.n=n;

      for(double U=0.5; U<4.0; U+=0.5)
      { char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
 
        iachm.SetParams(U, T, 0.5);
  
  
        //--------- RUN THE CODE---------//
        iachm.Run(&iaresult);


        // printout result
        char FN[300];
        sprintf(FN, "IACHM%s.IPT",bareFN);
        iaresult.PrintResult(FN);
        
        // use Pade to analytically continue Green's funtion to the real axis
        sprintf(FN, "Gw%s.IPT",bareFN);
        PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );
      }

    }
  }

}


