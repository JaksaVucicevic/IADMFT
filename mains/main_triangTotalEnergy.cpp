#include <limits>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"
#include "../source/IBZ.h"

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double U = atof(argv[1]);

  char EUFN[300];
  sprintf(EUFN,"E.U%.3f",U);
  FILE* EUFile = fopen(EUFN,"w");
  fclose(EUFile);


  int Nx = 500;
  int Ny = (int)( (double)Nx/(sqrt(3.0)) ); //Nx/2 if whole IBZ is used.
  printf("Nx: %d, Ny: %d\n",Nx,Ny);

  IBZ ibz(IBZtypes::TriangularLattice, Nx, Ny);

  int N = pow(2,13);

  double old_E=std::numeric_limits<double>::quiet_NaN();
  for(double T=0.01; T<0.5; T+=0.01)
  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);
 
    char bareFN[300];
    sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",0.5,U,T); 
    char FN[300];
    sprintf(FN, "IACHM%s",bareFN);
  
    iaresult.ReadFromFile(FN);

    complex<double> E = iaresult.InternalEnergy(T, &ibz);

    EUFile = fopen(EUFN,"a");
    fprintf(EUFile,"%.15le %.15le %.15le %.15le\n",T, real(E), imag(E), (real(E)-old_E)/0.01);
    fclose(EUFile);

    old_E = real(E);
  }
}

