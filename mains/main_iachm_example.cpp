#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"
#include <cstdlib>
#include <complex>
#include <cmath>

int main()
{
  int N = pow(2,13);

  IACHM iachm;

  iachm.SetBroydenOptions(false, false, 0); 
  iachm.SetLoopOptions(100, 1e-9);
  
  double T = 0.1;
  double U = 0.5;
  double n = 0.5;
  IAGRID iagrid(N,N, T);
  IAResult iaresult(&iagrid);
  iaresult.n=n; 
  iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau);
  iachm.PHSymmetricCase = true;
      
  char bareFN[300];
  sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
  iachm.SetParams(U, T, 0.5);
  iachm.LC->SetDoOutput(false);         
  //iachm.LC->SetOutputFileName(bareFN);         
  
  iachm.Run(&iaresult);
  char FN[300];
  sprintf(FN, "IACHM%s",bareFN);
  iaresult.PrintResult(FN);

  int Nnu = 5;
  double* nu = new double[Nnu];
  complex<double>* sigma = new complex<double>[Nnu];
  
  for(int n=1; n<=Nnu; n++)
  { printf("n: %d\n",n);
    sigma[n-1] = iaresult.OpticalConductivity(n);
    nu[n-1] = 2.0*n*pi*T;
  }
  sprintf(FN, "sigma%s",bareFN);
  PrintFunc(FN, Nnu, sigma, nu);

  delete [] nu;
  delete [] sigma;
}
