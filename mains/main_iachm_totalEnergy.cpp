#include <limits>
#include <cstdio>
#include <cstdlib>
#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"

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
  iachm.SetLoopOptions(100, 1e-9);

  // iterate over parameters
  for(double U=0.5; U<1.9; U+=0.5)
  {
    char EUFN[300];
    sprintf(EUFN,"E.U%.3f",U);
    FILE* EUFile = fopen(EUFN,"w");
    fclose(EUFile);

    double old_E=std::numeric_limits<double>::quiet_NaN();
    double old_smiE=std::numeric_limits<double>::quiet_NaN();  

    double Tstep=0.01;
    for(double T=0.01; T<0.801; T+=Tstep)
    {
      // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
      IAGRID iagrid(N,N, T);
      // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
      IAResult iaresult(&iagrid);
      // initialize FFT. this is also dependent on temperature
      iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); 
   
      iachm.PHSymmetricCase = true;

      iaresult.n=0.5;  

      iachm.SetParams(U, T, 0.5);
      
      //--------- RUN THE CODE---------//
      iachm.Run(&iaresult);
       
      complex<double> E = iaresult.InternalEnergy(T);
      complex<double> smiE = iaresult.SmartInternalEnergy(T);

      EUFile = fopen(EUFN,"a");
      fprintf(EUFile,"%.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", T, real(E), imag(E), real(smiE), imag(smiE), (real(E)-old_E)/Tstep, (real(smiE)-old_smiE)/Tstep);
      fclose(EUFile);

      old_E = real(E);
      old_smiE = real(smiE);

      // printout result
      char bareFN[300];
      sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",0.5,U,T);
      char FN[300];
      sprintf(FN, "IACHM%s",bareFN);
      iaresult.PrintResult(FN);
     
    }
  }

}

