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
  double T = atof(argv[1]);
  int Nx = 500;
  int Ny = (int)( (double)Nx/(sqrt(3.0)) ); //Nx/2 if whole IBZ is used.
  printf("Nx: %d, Ny: %d\n",Nx,Ny);

  IBZ ibz(IBZtypes::TriangularLattice, Nx, Ny);

  int N = pow(2,13);

  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);

    char rhoTFN[300];
    sprintf(rhoTFN,"rho.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);

    for(double U=1.0; U<5.0; U+=0.1)
    {   
      char bareFN[300];
      sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",0.5,U,T); 
      char FN[300];
      sprintf(FN, "IACHM%s",bareFN);
  
      iaresult.ReadFromFile(FN);

//      char rrFN[300];
//      sprintf(rrFN,"resRead%s",bareFN);
//      iaresult.PrintResult(rrFN);  
      int Nnu = 1000;
      double* nu = new double[Nnu];
      //complex<double>* Lambda = new complex<double>[Nnu];
      complex<double>* sigma = new complex<double>[Nnu];
      complex<double>* rho = new complex<double>[Nnu];
      complex<double>* logrho = new complex<double>[Nnu];

      for(int n=1; n<=Nnu; n++)
      { printf("n: %d\n",n);
        //Lambda[n] = iaresult.Lambda(n,&ibz); 
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

      rhoTFile = fopen(rhoTFN,"a");
      fprintf(rhoTFile,"%.15le %.15le %.15le %.15le\n", U, 1.0/real(sigma_w[0]), real(rho_w[0]), pow(10.0,real(logrho_w[0])));
      fclose(rhoTFile);

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
}

