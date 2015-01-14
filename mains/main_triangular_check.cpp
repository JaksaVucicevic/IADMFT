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
  if(argc<3) exit(0);
  double T = atof(argv[2]);
  double U = atof(argv[1]);

  int Nw = 1000;
  double* NIDOS = new double[Nw];
  double* w = new double[Nw];
  ReadFunc("triangular_dos", Nw, NIDOS, w);

  int Nx = 100;
  int Ny = (int)( (double)Nx/(sqrt(3.0)) ); //Nx/2 if whole IBZ is used.
  printf("Nx: %d, Ny: %d\n",Nx,Ny);

  IBZ ibz(IBZtypes::TriangularLattice, Nx, Ny);

  // initialize the object that does the calculation
  IACHM iachm;

  iachm.SetUseBethe(false);
  iachm.SetNIDOS(Nw, NIDOS, w);

  iachm.PHSymmetricCase = false;
  // set DMFT loop option
  iachm.SetBroydenOptions(false, false, 0); 
  iachm.SetLoopOptions(100, 1e-9);

  iachm.SetParams(U, T, -1000);

  int N = pow(2,13);

  IAGRID iagrid(N,N, T);
  IAResult iaresult(&iagrid);
  iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); 
  iaresult.n=0.5;
  iachm.Run(&iaresult);
  printf("iaresult.mu: %.3f\n",iaresult.mu);

  char bareFN[300];
  sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",0.5,U,T); 
  char FN[300];
  sprintf(FN, "IACHM%s",bareFN);
 
  iaresult.PrintResult(FN);
        
  // use Pade to analytically continue Green's funtion to the real axis
  sprintf(FN, "Gw%s",bareFN);
  PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );

  sprintf(FN, "Sigw%s",bareFN);
  PadeToFile(2000, iaresult.Sigma,  iaresult.omega, FN, 600, 4.0 );

  sprintf(FN, "Deltaw%s",bareFN);
  PadeToFile(2000, iaresult.Delta,  iaresult.omega, FN, 600, 4.0 );
if (false)
{
  int Nnu = 100;
  double* nu = new double[Nnu];
  complex<double>* Lambda = new complex<double>[Nnu+1];
  complex<double>* sigma = new complex<double>[Nnu];
  complex<double>* rho = new complex<double>[Nnu];
  complex<double>* logrho = new complex<double>[Nnu];

  Lambda[0] = iaresult.Lambda(0,&ibz);
  for(int n=1; n<=Nnu; n++)
  { printf("n: %d\n",n);
    Lambda[n] = iaresult.Lambda(n,&ibz); 
    sigma[n-1] = iaresult.OpticalConductivity(n, &ibz);
    rho[n-1] = 1.0/sigma[n-1];
    logrho[n-1] = log10(rho[n-1]);
    nu[n-1] = 2.0*n*pi*T;
  }
  sprintf(FN, "sigma%s",bareFN);
  PrintFunc(FN, Nnu, sigma, nu);
      
  int Nw = 200;
  double wmax =2.0;
  complex<double>* sigma_w = new complex<double>[Nw];
  complex<double>* rho_w = new complex<double>[Nw];
  complex<double>* logrho_w = new complex<double>[Nw];
  double* w = new double[Nw];
  for(int i=0; i<Nw; i++)
    w[i] = i * wmax/(Nw-1.0); 
  for(int n=0; n<Nnu; n++) sigma[n] = complex<double>(real(sigma[n]), 0);
  pade( Nnu, nu, sigma, 
        Nw,  w,  sigma_w );
  pade( Nnu, nu, rho, 
        Nw,  w,  rho_w );
  pade( Nnu, nu, logrho, 
        Nw,  w,  logrho_w );
  sprintf(FN, "sigma_w%s",bareFN);
  PrintFunc(FN, Nw, sigma_w, w);
  sprintf(FN, "rho_w%s",bareFN);
  PrintFunc(FN, Nw, rho_w, w);
  sprintf(FN, "logrho_w%s",bareFN);
  PrintFunc(FN, Nw, logrho_w, w);

  delete [] nu;
  delete [] sigma;
  delete [] rho;
  delete [] logrho;

  delete [] w; 
  delete [] sigma_w;
  delete [] rho_w;
  delete [] logrho_w;
}
}

