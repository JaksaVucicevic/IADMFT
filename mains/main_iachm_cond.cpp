#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/routines.h"

#include <cstdlib>
#include <cstdio>

double v(double eps) //this is actually v squared
{
  return 1.0-sqr(eps);
}

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double T = atof(argv[1]);
  // number of matsubara frequences to be used.
  // the fast fourier transform requires that this is a power of 2,
  // and that the number of tau points is the same
  int N = pow(2,13);

  
  // iterate over parameters
  //for(double T=0.01; T<0.51; T*=2.0)
  {
    // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
    IAGRID iagrid(N,N, T);
    // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
    IAResult iaresult(&iagrid);
    // initialize FFT. this is also dependent on temperature
   

    for(double n=0.5; n<0.7; n+=0.02)
    {
      char rhoTFN[300];
      sprintf(rhoTFN,"rho.n%.3f.T%.3f",n,T);
      FILE* rhoTFile = fopen(rhoTFN,"w");
      fclose(rhoTFile);

      for(double U=0.5; U<4.6; U+=0.5)
      { char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);

        // read result
        char FN[300];
        sprintf(FN, "IACHMs/IACHM%s",bareFN);

        iaresult.ReadFromFile(FN);

        // calculate sigma
        int Nnu = 300;
        double* nu = new double[Nnu];
        complex<double>* sigma = new complex<double>[Nnu];
        for(int m=1; m<=Nnu; m++)
        {  sigma[m-1] = iaresult.OpticalConductivity(m, &v);
           nu[m-1] = 2.0*m*pi*T;
        }
        sprintf(FN, "sigma%s",bareFN);
        PrintFunc(FN, Nnu, sigma, nu);        

        // pade sigma
        int Nw = 300;
        double wmax = 4.0;
        complex<double>* sigma_w = new complex<double>[Nw];
        double* w = new double[Nw];
        for(int i=0; i<Nw; i++)
          w[i] = i * wmax/(Nw-1.0); 
        pade( Nnu, nu, sigma, 
              Nw,  w,  sigma_w );
      
        sprintf(FN, "sigma_w%s",bareFN);
        PrintFunc(FN, Nw, sigma_w, w);

        rhoTFile = fopen(rhoTFN,"a");
        fprintf(rhoTFile,"%.15le %.15le\n", U, 1.0/real(sigma_w[0]));
        fclose(rhoTFile);

        delete [] nu;
        delete [] sigma;
      
        delete [] w; 
        delete [] sigma_w;
        
      }

    }
  }

}
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double T = atof(argv[1]);
  int Nx = 200;
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
      sprintf(FN, "IACHMs/IACHM%s",bareFN);
  
      iaresult.ReadFromFile(FN);

//      char rrFN[300];
//      sprintf(rrFN,"resRead%s",bareFN);
//      iaresult.PrintResult(rrFN);  
      int Nnu = 100;
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
*/
