#include <cstdlib>
#include <complex>
#include <cmath>
#include <cstdio>
#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"
#include "../source/IBZ.h"


double TriangularConductivity(int Nw, complex<double>* Sigma, double* w, double mu,
                              double T, int Nkx, int Nky, int Nnu, double numax)
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  IBZ ibz(IBZtypes::TriangularLattice, Nkx, Nky );
  printf("-----ibz ready\n");
  complex<double> sum = 0.0;
  
  double k = 20.0;		//determines the nu range

  double nu_start = -k*T;
  double nu_end = k*T;

  if (nu_start < - numax) nu_start = - numax;
  if (nu_end > numax) nu_end =  numax;


  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { printf("-- nu: %.3f\n",nu);

    complex<double> Sigma_nu=interpl(Nw, Sigma, w, nu);
    double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                              ) 
                                        )   
                           );

    for(int i=0; i<Nkx; i++)
    for(int j=0; j<Nky; j++)
    { 
      double v = ibz.velocity[i][j];
 
      complex<double> G = 1.0/( nu + mu - ibz.epsilon[i][j] - Sigma_nu);
      double rho = - imag(G) / pi;

      ibz.summand[i][j] = sqr(rho) * sqr(v) * fprim;

      //if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    sum+=ibz.sum();

    //if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  //if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * real(sum) * dnu ;
}



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

  /*  char rhoTFN[300];
    sprintf(rhoTFN,"RArhos/rho.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);
*/
    for(double U=0.5; U<5.0; U+=10000.1)
    {   
      char bareFN[300];
      sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",0.5,U,T); 

      //--------- read data -----------//
      char FN[300];
      sprintf(FN, "IACHMs/IACHM%s",bareFN);
  
      iaresult.ReadFromFile(FN);
      printf(" --- read: %s\n", FN);  
     
      //--------- pade self-energy -----------//
      int Nw = 1000;
      double wmax =5.0;
      complex<double>* Sigma_w = new complex<double>[Nw];
      double* w = new double[Nw];
      for(int i=0; i<Nw; i++)
        w[i] = - wmax+ 2.0 * i * wmax/(Nw-1.0); 

      pade( 2000, iaresult.omega, iaresult.Sigma, 
            Nw,  w,  Sigma_w );

      printf("---------------- did pade on Sigma\n");

      char SFN[300];
      sprintf(SFN, "Sigws/Sigw%s",bareFN);
      PrintFunc(SFN, Nw, Sigma_w, w);

      printf("---------------- printed to file\n");

      //------- pade green's function ---------//
      complex<double>* G_w = new complex<double>[Nw];

      pade( 2000, iaresult.omega, iaresult.G, 
            Nw,  w,  G_w );

      printf("---------------- did pade on G\n");

      char GFN[300];
      sprintf(GFN, "Gws/Gw%s",bareFN);
      PrintFunc(GFN, Nw, G_w, w);

      printf("---------------- printed to file\n");
      
      //-------- conductivity from real axis data -----------//      
/*
      double sigma = TriangularConductivity(Nw, Sigma_w, w, iaresult.mu,
                                            T, Nx, Ny, 400, 10000);

      rhoTFile = fopen(rhoTFN,"a");
      fprintf(rhoTFile,"%.15le %.15le\n", U, 1.0/sigma);
      fclose(rhoTFile);

      printf("------- resistivity done!\n");
*/

      delete [] w; 
      delete [] G_w;
      delete [] Sigma_w;
    }   
  }
}

