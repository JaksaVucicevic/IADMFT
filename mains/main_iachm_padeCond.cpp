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

double phi_bethe(double eps)
{
  return pow( 1.0-sqr(eps),  3.0 );
}


double Conductivity(int N, complex<double>* Sigma, double* omega, double T, double mu, int Neps, int Nnu, const char * integrandFN, bool excludeSmallOmega, double (*phi)(double))
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");

  double sum = 0.0;
  
  double k = 20.0;		//determines the nu range
  double W = 10.0;		//determines the epsilon range
  double NIDOS_EDGE = 1.0;	//edge of the non-interacting band

  //double dnu = 2.0 * k * T / (double) Nnu;
  //for (double nu = -k*T; nu < k*T; nu += dnu )
 
  double nu_start = -k*T;
  double nu_end = k*T;

  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { 
    //---------- ONLY FOR INSULATOR -----------//
    if ((abs(nu)<0.1)and(excludeSmallOmega))
      continue;
    //-----------------------------------------//

    complex<double> Sigma_nu=interpl(N, Sigma, omega, nu);
    double eps_center = nu+mu-real(Sigma_nu);
    double eps_start=eps_center-W*abs(imag(Sigma_nu)); 
    double eps_end=eps_center+W*abs(imag(Sigma_nu));
    if (eps_start < -NIDOS_EDGE) eps_start = -1.0;
    if (eps_end > NIDOS_EDGE)    eps_end = 1.0;
    double deps = (eps_end-eps_start) / Neps ;
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      complex<double> G = 1.0/( nu + mu - eps - Sigma_nu);
      double rho = - imag(G) / pi;
      double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                                ) 
                                          )   
                           );

      double integrand = sqr(rho) * (*phi)(eps) * fprim; // only for bethe lattice! v = rho
      sum += integrand * deps;
      if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, (*phi)(eps));  
    }
    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}



int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double T = atof(argv[1]);

  int N = pow(2,13);

  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);

    char cmd[300];
    sprintf(cmd,"mkdir paderhos; mkdir Sigws; mkdir Gws");
    system(cmd);

    for(double n=0.5; n<0.7; n+=0.02)
    {
      char rhoTFN[300];
      sprintf(rhoTFN,"paderhos/rho.n%.3f.T%.3f",n,T);
      FILE* rhoTFile = fopen(rhoTFN,"w");
      fclose(rhoTFile);

      for(double U=0.5; U<4.6; U+=0.5)
      {   
        char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T); 

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

        double sigma =  Conductivity(Nw, Sigma_w, w, T, iaresult.mu, 400, 400, NULL, false, &phi_bethe);

        rhoTFile = fopen(rhoTFN,"a");
        fprintf(rhoTFile,"%.15le %.15le\n", U, 1.0/sigma);
        fclose(rhoTFile);

        printf("------- resistivity done!\n");

        delete [] w; 
        delete [] G_w;
        delete [] Sigma_w;
      } 
    }  
  }
}

