#include "../source/IAresArray.h"
#include "../source/IAResult.h"
#include "../source/IAGRID.h"
#include "../source/arrayInitializers.h"
#include "../source/ctqmcDHM.h"
#include "../source/ctqmcSIAM.h"
#include "../source/routines.h"
#include "../source/pade.h"
#include <omp.h>

#include <cstdlib>
#include <cstdio>

double t = 1.0/6.0;

double epsilon(double kx, double ky, double kz)
{
  return -2.0*t*(cos(kx)+cos(ky)+cos(kz));
}

double velocity(double kx)
{
  return 2.0*t*sin(kx);
}

//current current correlation
complex<double> Lambda(int N, complex<double>* Sigma, double* omega, double mu, int Nk, int n)
{
  //-- bosonic frequency
 
  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma[i];
    Sig_large[N-1-i] = conj(Sigma[i]); // CONJ: G(iw)=G*(-iw) 
  }

  complex<double> sum = 0.0;
  int Nsum = 200;
  //for (int m=0; m<2*N-n; m++)
  for (int m=N-1-Nsum; m<N+Nsum; m++)
  { 
    complex<double> ksum=0.0;
    for(int i=0; i<Nk; i++)
    for(int j=0; j<Nk; j++)
    for(int l=0; l<Nk; l++)
    { 
      double kx = i*pi/(double)(Nk);
      double ky = j*pi/(double)(Nk);
      double kz = l*pi/(double)(Nk);

      double v = velocity(kx);
      double eps = epsilon(kx,ky,kz);

      ksum  +=  sqr(v) 
              * 2.0 
              * ( 1.0
                  / ( ( ii*iw_large[m] + mu - eps - Sig_large[m] ) 
                       * 
                      ( ii*iw_large[m+n] + mu - eps - Sig_large[m+n] ) 
                    )
                );
    }
    sum += ksum/(double)pow(Nk,3);
  }

  delete [] iw_large;
  delete [] Sig_large;

  return sum;
}

complex<double> OpticalConductivity(int N, complex<double>* Sigma, double* omega, double mu, int Nk, int n, complex<double> Lambda0)
{
  return (1.0/(2.0*pi*n))* ( Lambda(N, Sigma, omega, mu, Nk, n) - Lambda0 );
}


double CubicConductivity(int Nw, complex<double>* Sigma, double* w, double mu,
                         double T, int Nk, int Nnu, double numax, bool excludeSmallOmega=false)
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  complex<double> sum = 0.0;
  
  double k = 20.0;		//determines the nu range

  double nu_start = -k*T;
  double nu_end = k*T;

  if (nu_start < - numax) nu_start = - numax;
  if (nu_end > numax) nu_end =  numax;


  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { printf("-- nu: %.3f\n",nu);

    //---------- ONLY FOR INSULATOR -----------//
    if ((abs(nu)<0.1)and(excludeSmallOmega))
      continue;
    //-----------------------------------------//

    complex<double> Sigma_nu=interpl(Nw, Sigma, w, nu);
    double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                              ) 
                                        )   
                           );

    double ksum=0.0;
    for(int i=0; i<Nk; i++)
    for(int j=0; j<Nk; j++)
    for(int l=0; l<Nk; l++)
    { 
      double kx = i*pi/(double)(Nk);
      double ky = j*pi/(double)(Nk);
      double kz = l*pi/(double)(Nk);

      double v = velocity(kx);
      double eps = epsilon(kx,ky,kz);
 
      complex<double> G = 1.0/( nu + mu - eps - Sigma_nu);
      double rho = - imag(G) / pi;

      ksum += sqr(rho) * sqr(v) * fprim;

      //if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    sum += (ksum/(double)pow(Nk,3));

    //if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  //if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * real(sum) * dnu ;
}


//=================================== SMART HISTOGRAMS ============================================//

void ClearSmartHistogramOutput(const char* sufix, const char* sufix2, double U, double W, double T)
{
  if (W!=0.0)
  { char vsUFN[300];
    sprintf(vsUFN, "%s%s_smartHistogram_vs_U.W%.3f.T%.3f", sufix, sufix2, W, T);
    FILE* vsUFile = fopen(vsUFN,"w");
    fclose(vsUFile);   
  }

  if (U!=0.0)
  { char vsWFN[300];
    sprintf(vsWFN, "%s%s_smartHistogram_vs_W.U%.3f.T%.3f",sufix, sufix2, U, T);
    FILE* vsWFile = fopen(vsWFN,"w");
    fclose(vsWFile); 
  }
}

void AddNewLine(const char* sufix, const char* sufix2, double W, double T)
{
  char vsUFN[300];
  sprintf(vsUFN, "%s%s_smartHistogram_vs_U.W%.3f.T%.3f", sufix, sufix2, W, T);
  FILE* vsUFile = fopen(vsUFN,"a");
  fprintf(vsUFile,"\n");
  fclose(vsUFile);   
}

void ReadDataAndPrintSmartHistogram(const char* FN, const char* sufix, const char* sufix2, int N, double* Y, double* X, int Nbins, double* P, double* y, double U, double W, double T, bool Logarithmic=false)
{
  char fn[300];
  sprintf(fn,"%s/%ss%s",FN,sufix,sufix2);

  if (FileExists(fn)) ReadFunc(fn, N, Y, X); 
  else { printf(">>>>>>>>>> FILE NOT FOUND: %s\n\n",fn); return;  }

  SmartHistogram(N, Y, X, Nbins, y, P, Logarithmic);

  char fn2[300];
  sprintf(fn2,"%s/%s%s_smartHistogram",FN,sufix,sufix2);
  PrintFunc(fn2, Nbins, P, y);     
  printf("printed %s%s smartHistogram!\n",sufix,sufix2);

  char vsUFN[300];
  sprintf(vsUFN, "%s%s_smartHistogram_vs_U.W%.3f.T%.3f", sufix, sufix2, W, T);
  FILE* vsUFile = fopen(vsUFN,"a");
  for(int b=0; b<Nbins; b++)
    fprintf(vsUFile,"%.15le %.15le %.15le\n", U, y[b], P[b]);
  fclose(vsUFile);   

  char vsWFN[300];
  sprintf(vsWFN, "%s%s_smartHistogram_vs_W.U%.3f.T%.3f",sufix, sufix2, U, T);
  FILE* vsWFile = fopen(vsWFN,"a");
  for(int b=0; b<Nbins; b++)
    fprintf(vsWFile,"%.15le %.15le %.15le\n", U, y[b], P[b]);
  fprintf(vsWFile,"\n");
  fclose(vsWFile); 
}

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];

  //for(double T=0.1; T>0.0009; T/=10.0)
  {      
    for(double W = 0.25; W<12.1; W+=0.25)
    {
      ClearSmartHistogramOutput("Aw0", "", 0.0, W, T);
      ClearSmartHistogramOutput("Aw0", "_pade", 0.0, W, T);
      ClearSmartHistogramOutput("n", "", 0.0, W, T);
      ClearSmartHistogramOutput("tau", "", 0.0, W, T);
      ClearSmartHistogramOutput("tau", "_pade", 0.0, W, T);
      ClearSmartHistogramOutput("tau", "_RASC", 0.0, W, T);
      ClearSmartHistogramOutput("tau", "_RASC_pade", 0.0, W, T);
      ClearSmartHistogramOutput("rmu", "", 0.0, W, T);
      ClearSmartHistogramOutput("ReSigw0", "", 0.0, W, T);
      ClearSmartHistogramOutput("Z", "", 0.0, W, T);
      ClearSmartHistogramOutput("Z", "_RASC", 0.0, W, T);
    }

    for(double U = 0.25; U<9.6; U+=0.25)
    {
      ClearSmartHistogramOutput("Aw0", "", U, 0.0, T);
      ClearSmartHistogramOutput("Aw0", "_pade", U, 0.0, T);
      ClearSmartHistogramOutput("n", "", U, 0.0, T);
      ClearSmartHistogramOutput("tau", "", U, 0.0, T);
      ClearSmartHistogramOutput("tau", "_pade", U, 0.0, T);
      ClearSmartHistogramOutput("tau", "_RASC", U, 0.0, T);
      ClearSmartHistogramOutput("tau", "_RASC_pade", U, 0.0, T);
      ClearSmartHistogramOutput("rmu", "", U, 0.0, T);
      ClearSmartHistogramOutput("ReSigw0", "", U, 0.0, T);
      ClearSmartHistogramOutput("Z", "", U, 0.0, T);
      ClearSmartHistogramOutput("Z", "_RASC", U, 0.0, T);
    }

    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { 
      for(double W = 0.25; W<12.1; W+=0.25)
      { //==========//

        //load the result and clear previous sorted files
        sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);    
      
        printf("======================= WORKING: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);    
      
        //------- prepare histogram arrays ------//
        int Nbins = 80;
        double* P = new double[Nbins];
        double* y = new double[Nbins];

        double* X = new double[Nsites];     
        double* Y = new double[Nsites];

        //--------- Aw0s
        ReadDataAndPrintSmartHistogram(FN, "Aw0", "", Nsites, Y, X, Nbins, P, y, U, W, T);
      
        //--------- Aw0s_pade
        ReadDataAndPrintSmartHistogram(FN, "Aw0", "_pade", Nsites, Y, X, Nbins, P, y, U, W, T);

        //--------- taus
        ReadDataAndPrintSmartHistogram(FN, "tau", "", Nsites, Y, X, Nbins, P, y, U, W, T, true);

        //--------- taus_pade
        ReadDataAndPrintSmartHistogram(FN, "tau", "_pade", Nsites, Y, X, Nbins, P, y, U, W, T, true);

        //--------- taus_RASC
        ReadDataAndPrintSmartHistogram(FN, "tau", "_RASC", Nsites, Y, X, Nbins, P, y, U, W, T, true);

        //--------- taus_RASC_pade
        ReadDataAndPrintSmartHistogram(FN, "tau", "_RASC_pade", Nsites, Y, X, Nbins, P, y, U, W, T, true);

        //--------- n
        ReadDataAndPrintSmartHistogram(FN, "n", "", Nsites, Y, X, Nbins, P, y, U, W, T);  

        //--------- rmu
        ReadDataAndPrintSmartHistogram(FN, "rmu", "", Nsites, Y, X, Nbins, P, y, U, W, T);  

        //--------- ReSigw0
        ReadDataAndPrintSmartHistogram(FN, "ReSigw0", "", Nsites, Y, X, Nbins, P, y, U, W, T);  

        //--------- Z
        ReadDataAndPrintSmartHistogram(FN, "Z", "", Nsites, Y, X, Nbins, P, y, U, W, T);  

        //--------- Z_RASC
        ReadDataAndPrintSmartHistogram(FN, "Z", "_RASC", Nsites, Y, X, Nbins, P, y, U, W, T);  

        delete [] P;
        delete [] y;
        delete [] X;
        delete [] Y;
        
      }
      for(double W = 0.25; W<12.1; W+=0.25)
      {
        AddNewLine("Aw0", "",  W, T);
        AddNewLine("Aw0", "_pade",  W, T);
        AddNewLine("n", "",  W, T);
        AddNewLine("tau", "",  W, T);
        AddNewLine("tau", "_pade",  W, T);
        AddNewLine("tau", "_RASC",  W, T);
        AddNewLine("tau", "_RASC_pade",  W, T);
        AddNewLine("rmu", "",  W, T);
        AddNewLine("ReSigw0", "",  W, T);
        AddNewLine("Z", "", W, T);
        AddNewLine("Z", "_RASC",  W, T);
      }    
    }
  } 

  return 0;
}




//============================================== read Aw0s, calculate avg, typical and mpv, 3d and 2d printout UxW ==============================//
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    // prepare files
    
    char charAw0sTFN[300];
    sprintf(charAw0sTFN,"charAw0s.T%.3f",T);
    FILE* charAw0sTFile = fopen(charAw0sTFN,"w");
    fclose(charAw0sTFile); 
    
    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char charAw0sWFN[300];
      sprintf(charAw0sWFN, "charAw0s_vs_U.W%.3f.T%.3f",W,T);
      FILE* charAw0sWFile = fopen(charAw0sWFN,"w");
      fclose(charAw0sWFile);
    }
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    {

      char charAw0sUFN[300];
      sprintf(charAw0sUFN, "charAw0s_vs_W.U%.3f.T%.3f",U,T);
      FILE* charAw0sUFile = fopen(charAw0sUFN,"w");
      fclose(charAw0sUFile);

 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//
      double* Aw0s = new double[Nsites];
      double* epsilons = new double[Nsites];

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f/Aw0s",W,T,U);
      if (not FileExists(FN)) continue;
      ReadFunc(FN, Nsites, Aw0s, epsilons);     
     
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);
     
      double Aw0_typ = typical(Nsites, Aw0s);
      double Aw0_avg = average(Nsites, Aw0s);
      
      int Nbins = 12;
      double* x = new double[Nbins];
      double* P = new double[Nbins];
      Histogram(Nsites, Aw0s, Nbins, x, P);      

      int index = 0;
      maximum(Nbins, P, &index);
      double Aw0_mpv = x[index];

      delete [] Aw0s;
      delete [] epsilons;
      delete [] x;
      delete [] P;

      printf("----------> calculated characteristic values, now about to printout results");

      charAw0sUFile = fopen(charAw0sUFN,"a");
      fprintf(charAw0sUFile,"%.15le %.15le %.15le %.15le\n", W, Aw0_avg, Aw0_typ, Aw0_mpv);
      fclose(charAw0sUFile);

      charAw0sTFile = fopen(charAw0sTFN,"a");
      fprintf(charAw0sTFile,"%.15le %.15le %.15le %.15le %.15le\n", U, W, Aw0_avg, Aw0_typ, Aw0_mpv);
      fclose(charAw0sTFile);

      char charAw0sWFN[300];
      sprintf(charAw0sWFN, "charAw0s_vs_U.W%.3f.T%.3f",W,T);
      FILE* charAw0sWFile = fopen(charAw0sWFN,"a");
      fprintf(charAw0sWFile,"%.15le %.15le %.15le %.15le\n", U, Aw0_avg, Aw0_typ, Aw0_mpv);
      fclose(charAw0sWFile);
      
    }
    charAw0sTFile = fopen(charAw0sTFN,"a");
    fprintf(charAw0sTFile,"\n");
    fclose(charAw0sTFile);

    }
  }

  return 0;
}

*/

//=================================== ALL vs T ============================================//
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double W = atof(argv[1]);
//  double T = atof(argv[1]);
//  double U = atof(argv[2]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { 
    // prepare files
    
//    char rhoWFN[300];
//    sprintf(rhoWFN,"rho.W%.3f",W);
//    FILE* rhoWFile = fopen(rhoWFN,"w");
//    fclose(rhoWFile);

//    char rhoRASCWFN[300];
//    sprintf(rhoRASCWFN,"rho_RASC.W%.3f",W);
//    FILE* rhoRASCWFile = fopen(rhoRASCWFN,"w");
//    fclose(rhoRASCWFile);
       
    // ===================== ITERATE ======================== //
    //for(double U = 1.9; U<W; U+=0.5)
    for(double U = W-1.0; U<W+1.0; U+=0.05)
    { //prepare files for which names need to be prepared only once
      char Aw0hFN[300];
      sprintf(Aw0hFN, "Aw0_histogram_vs_T.W%.3f.U%.3f",W,U);
      FILE* Aw0hFile = fopen(Aw0hFN,"w");
      fclose(Aw0hFile);

      char nhFN[300];
      sprintf(nhFN, "n_histogram_vs_T.W%.3f.U%.3f",W,U);
      FILE* nhFile = fopen(nhFN,"w");
      fclose(nhFile); 

      char tauhFN[300];
      sprintf(tauhFN, "tau_histogram_vs_T.W%.3f.U%.3f",W,U);
      FILE* tauhFile = fopen(tauhFN,"w");
      fclose(tauhFile);

      char rmuhFN[300];
      sprintf(rmuhFN, "rmu_histogram_vs_T.W%.3f.U%.3f",W,U);
      FILE* rmuhFile = fopen(rmuhFN,"w");
      fclose(rmuhFile);

      char ReSigw0hFN[300];
      sprintf(ReSigw0hFN, "ReSigw0_histogram_vs_T.W%.3f.U%.3f",W,U);
      FILE* ReSigw0hFile = fopen(ReSigw0hFN,"w");
      fclose(ReSigw0hFile);

      char rhoWUFN[300];
      sprintf(rhoWUFN,"rho.W%.3f.U%.3f",W,U);
      FILE* rhoWUFile = fopen(rhoWUFN,"w");
      fclose(rhoWUFile);

      char rhoRASCWUFN[300];
      sprintf(rhoRASCWUFN,"rho_RASC.W%.3f.U%.3f",W,U);
      FILE* rhoRASCWUFile = fopen(rhoRASCWUFN,"w");
      fclose(rhoRASCWUFile);

 
    for(double T = 0.002; T<0.55; T*=2.0)
    { //==========//

      IAGRID g( Nw, Nw, T );

      //all the storage arrays
      IAresArray a( Nsites, &g );

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      char cmd[300];
      //sprintf(cmd,"rm %s/mu*",FN);
      //system(cmd);
      
      a.ReadFromUberMinimalFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);
     
      //prints ns and ns_sorted
      a.Print_ns_and_mus(FN);

      double* ns = a.Get_ns();     
      sprintf(fn,"%s/n_histogram",FN);
      a.Print_Histogram(fn, ns);
      printf("printed n histogram!\n");

      printf("printed ns and mus!\n");

      
      //--- SPLINE

      double* Aw0s = a.Get_Aw0s();
      double* taus = a.Get_taus();  
      double* mus = a.Get_mus();
      double* rmus = a.Get_renormalized_mus();
      double* ReSigw0s = a.Get_ReSigw0s();

      //print sorted
      sprintf(fn,"%s/Aw0s",FN);
      PrintFunc(fn, Nsites, Aw0s, mus);

      printf("printed Aw0s!\n");

      sprintf(fn,"%s/taus",FN);
      PrintFunc(fn, Nsites, taus, mus);
      printf("printed taus!\n");
      
      //print histograms
      sprintf(fn,"%s/Aw0_histogram",FN);
      a.Print_Histogram(fn, Aw0s);
      printf("printed Aw0 histogram!\n");

      sprintf(fn,"%s/tau_histogram",FN);
      a.Print_Histogram(fn, taus, true);
      printf("printed tau histogram!\n");

      sprintf(fn,"%s/rmu_histogram",FN);
      a.Print_Histogram(fn, rmus);
      printf("printed renormalized mu histogram!\n");

      sprintf(fn,"%s/ReSigw0_histogram",FN);
      a.Print_Histogram(fn, ReSigw0s, true);
      printf("printed ReSigw0 histogram!\n");


      //--- PADE

      double* pade_Aw0s = a.Get_pade_Aw0s();
      double* pade_taus = a.Get_pade_taus();  

      //print sorted
      sprintf(fn,"%s/Aw0s_pade",FN);
      PrintFunc(fn, Nsites, pade_Aw0s, mus);
      printf("printed Aw0s!\n");

      sprintf(fn,"%s/taus_pade",FN);
      PrintFunc(fn, Nsites, pade_taus, mus);
      printf("printed taus!\n");
      
      //print histograms
      sprintf(fn,"%s/Aw0_pade_histogram",FN);
      a.Print_Histogram(fn, pade_Aw0s);
      printf("printed Aw0 histogram!\n");

      sprintf(fn,"%s/tau_pade_histogram",FN);
      a.Print_Histogram(fn, pade_taus, true);
      printf("printed tau histogram!\n");


     
      //------------add to the 3D histogram files--------//
      double* x;
      double* P;
      double* pade_x;
      double* pade_P;

      int Nbins = 12;

      //---Aw0s, both spline and pade
      a.Get_Histogram(Aw0s, Nbins, x, P);
      a.Get_Histogram(pade_Aw0s, Nbins, pade_x, pade_P);

      Aw0hFile = fopen(Aw0hFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(Aw0hFile,"%.15le %.15le %.15le %.15le %.15le\n",W,x[b],P[b],pade_x[b],pade_P[b]);
      fprintf(Aw0hFile,"\n");
      fclose(Aw0hFile);   

      delete [] x;
      delete [] P;
      delete [] pade_x;
      delete [] pade_P;
      printf("added to Aw0_histogram_vs_W and U!\n");

      //---taus, both spline and pade
      a.Get_Histogram(taus, Nbins, x, P, true);
      a.Get_Histogram(pade_taus, Nbins, pade_x, pade_P, true);

      tauhFile = fopen(tauhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhFile,"%.15le %.15le %.15le %.15le %.15le\n",W,x[b],P[b],pade_x[b],pade_P[b]);
      fprintf(tauhFile,"\n");
      fclose(tauhFile);   

      delete [] x;
      delete [] P;
      delete [] pade_x;
      delete [] pade_P;

      printf("added to tau_histogram_vs_W!\n");

      //---ns
      a.Get_Histogram(ns, Nbins, x, P);

      nhFile = fopen(nhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(nhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(nhFile,"\n");
      fclose(nhFile);   

      delete [] x;
      delete [] P;
      printf("added to n_histogram_vs_W!\n");

      //---Renomralized mus
      a.Get_Histogram(rmus, Nbins, x, P);

      rmuhFile = fopen(rmuhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(rmuhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(rmuhFile,"\n");
      fclose(rmuhFile);   

      delete [] x;
      delete [] P;
      printf("added to rmu_histogram_vs_W!\n");

      //---ReSigw0s
      a.Get_Histogram(ReSigw0s, Nbins, x, P);

      ReSigw0hFile = fopen(ReSigw0hFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ReSigw0hFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(ReSigw0hFile,"\n");
      fclose(ReSigw0hFile);   

      delete [] x;
      delete [] P;
      printf("added to rmu_histogram_vs_W!\n");

      delete [] ns;
      delete [] mus;
      delete [] Aw0s;
      delete [] taus;    
      delete [] pade_Aw0s;
      delete [] pade_taus;    

      //---------------- Gavg, Delta and Sigma effective----------------------//
      complex<double>* Sigma;
      complex<double>* G;
      complex<double>* Delta;
      double* iw;
     
      int M;
      sprintf(fn,"%s/Sigma_eff",FN);
      if (not ReadFunc(fn, M, Sigma, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) Sigma[i] = complex<double>(U/2.0,imag(Sigma[i]));
      //sprintf(fn,"%s/Sigma_eff_read",FN);
      //PrintFunc(fn, M, Sigma, iw);
      printf("------- read Sigma effective!\n\n");
      
      sprintf(fn,"%s/G_avg",FN);
      if (not ReadFunc(fn, M, G, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) G[i] = complex<double>(0.0,imag(G[i]));
      //sprintf(fn,"%s/G_avg_read",FN);
      //PrintFunc(fn, M, G, iw);
      printf("------- read G average!\n\n");

      sprintf(fn,"%s/Delta",FN);
      if (not ReadFunc(fn, M, Delta, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) Delta[i] = complex<double>(0.0,imag(Delta[i]));
      printf("------- read Delta!\n\n");

      int Mw = 2000;
      complex<double>* Gw = new complex<double>[Mw];
      complex<double>* Sigmaw = new complex<double>[Mw];
      complex<double>* Deltaw = new complex<double>[Mw];

      double* w = new double[Mw];
      double wmax = max(U/2.0,W/2.0) + 2.0;
      for(int i=0; i<Mw; i++)
        w[i] = - wmax + i * 2.0*wmax/(Mw-1.0);     

      //PrintFunc("w_created",Mw,w); 
      printf("------- about to do pade...\n\n");

     //this is already done, we don't need it
      int Nf = a.GetActualN();
      if (Nf>50) Nf = 50;

      pade( Nf, iw, G, 
            Mw, w, Gw );
      sprintf(fn,"%s/G_avg_w",FN);
      PrintFunc(fn, Mw, Gw, w);
   
      pade( Nf, iw, Sigma, 
            Mw, w, Sigmaw );
      sprintf(fn,"%s/Sigma_eff_w",FN);
      PrintFunc(fn, Mw, Sigmaw, w);

      pade( Nf, iw, Delta, 
            Mw, w, Deltaw );
      sprintf(fn,"%s/Delta_w",FN);
      PrintFunc(fn, Mw, Deltaw, w);

      printf(">>>>>>>> Did pade on G_avg, Delta and Sigma_eff and printed to files!\n\n");

      complex<double>* Sigma_eff_w_RASC = new complex<double>[Mw];
      for(int i=0; i<Mw; i++)
        Sigma_eff_w_RASC[i] = w[i]+U/2.0-Deltaw[i]-1.0/Gw[i];

      sprintf(fn,"%s/Sigma_eff_w_RASC",FN);
      PrintFunc(fn, Mw, Sigma_eff_w_RASC, w);

      printf(">>>>>>>> Printed Sigma_eff_w obtained from self-consistency on the real axis!\n\n");


      //-------- conductivity from imag axis Sigma effective -----------//      
      printf("------- sigma(nu)...\n");
      complex<double> Lambda0 = Lambda(M, Sigma, iw, U/2.0, 16, 0); //calculate zero frequency current-current correlation only once

      int Nn = 400;
      complex<double>* sigma = new complex<double>[Nn];
      double* nu = new double[Nn]; 
      for(int n=0; n<Nn; n++)
      {  sigma[n] = OpticalConductivity(M, Sigma, iw, U/2.0, 16, n+1, Lambda0);
         nu[n] = 2.0*(n+1)*pi*T;
      }
      
      sprintf(fn,"%s/sigma_nu",FN);
      PrintFunc(fn,Nn,sigma,nu);

      double rho = 1.0/real(pade( Nn, nu, sigma, 0.0 ));
      printf("------- sigma(nu) DONE!\n");

      //-------- conductivity from real axis Sigma effective -----------//      

      printf("------- resistivity...\n");
      double sgm =  CubicConductivity(Mw, Sigmaw, w, U/2.0, T, 16, 400, 1000.0);

//      rhoWFile = fopen(rhoWFN,"a");
//      fprintf(rhoWFile,"%.15le %.15le %.15le %.15le\n", U, T, 1.0/sgm, rho);
//      fclose(rhoWFile);

      rhoWUFile = fopen(rhoWUFN,"a");
      fprintf(rhoWUFile,"%.15le %.15le %.15le\n", T, 1.0/sgm, rho);
      fclose(rhoWUFile);

      //-------- conductivity from RASC Sigma effective -----------//       
      {
      printf("------- resistivity RASC...\n");
      double sgm =  CubicConductivity(Mw, Sigma_eff_w_RASC, w, U/2.0, T, 16, 400, 1000.0);
      double sgm2 =  CubicConductivity(Mw, Sigma_eff_w_RASC, w, U/2.0, T, 16, 400, 1000.0, true);

//      rhoRASCWFile = fopen(rhoRASCWFN,"a");
//      fprintf(rhoRASCWFile,"%.15le %.15le %.15le %.15le\n", U, T, 1.0/sgm, 1.0/sgm2);
//      fclose(rhoRASCWFile);

      rhoRASCWUFile = fopen(rhoRASCWUFN,"a");
      fprintf(rhoRASCWUFile,"%.15le %.15le %.15le\n", T, 1.0/sgm, 1.0/sgm2);
      fclose(rhoRASCWUFile);

      printf("------- resistivity done!\n");
      }

      delete [] G;
      delete [] Sigma;
      delete [] Gw;
      delete [] Sigmaw;
      delete [] iw;
      delete [] w;
      delete [] nu;

      //release H0     
    }
//    rhoWFile = fopen(rhoWFN,"a");
//    fprintf(rhoWFile,"\n");
//    fclose(rhoWFile);

//    rhoRASCWFile = fopen(rhoRASCWFN,"a");
//    fprintf(rhoRASCWFile,"\n");
//    fclose(rhoRASCWFile);

  }
  }

  return 0;
}
*/





//============================================== G_epsilon pade ================================================//
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );  

    
    // ===================== ITERATE ======================== //
    for(double U = 1.0; U<6.0; U+=1.0)
    { 
 
    for(double W = 1.0; W<9.0; W+=1.0)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      a.ReadFromUberMinimalFiles(FN);


      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);


      //Pade on Green's function and Self-energy
      for(int id=9; id<Nsites; id+=10)
      { char pFN[300];
        sprintf(pFN,"%s/Gw.%d",FN,id);     

        int M = a.GetActualN();
        if (M>1000) M=1000;
        PadeToFile( M, a.r[id].G,  a.r[id].omega, pFN, 500, U/2.0+max(2.0, W/2.0) );
        //sprintf(pFN,"%sSigw.%d",FN,id);
        //PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 500, 3.0 );

        printf("pade, site id=%d : DONE!\n",id);
      }

    }


  }
  }

  return 0;
}
*/


//============================================== Zs and Z histograms ================================================//
/*int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  char sufix[300];
  sprintf(sufix,"_RASC");

  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );  

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char ZhWFN[300];
      sprintf(ZhWFN, "Z%s_histogram_vs_U.W%.3f.T%.3f",sufix,W,T);
      FILE* ZhWFile = fopen(ZhWFN,"w");
      fclose(ZhWFile);
    }
    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { //prepare files for which names need to be prepared only once
      char ZhFN[300];
      sprintf(ZhFN, "Z%s_histogram_vs_W.U%.3f.T%.3f",sufix,U,T);
      FILE* ZhFile = fopen(ZhFN,"w");
      fclose(ZhFile);
 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      a.ReadFromUberMinimalFiles(FN);
      double* mus = a.Get_mus();

      complex<double>* Delta;
      double* iw;
      int Miw;
      sprintf(fn,"%s/Delta",FN);
      if (not ReadFunc(fn, Miw, Delta, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<Miw; i++) Delta[i] = complex<double>(0.0,imag(Delta[i]));

      printf("------- read Delta!\n\n");

      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);

      //double* Zs = a.Get_Zs();
      double* Zs = a.Get_Zs_RASC(Delta);

      sprintf(fn,"%s/Zs%s",FN,sufix);
      PrintFunc(fn, Nsites, Zs, mus);
      printf("printed Zs!\n");
      
      sprintf(fn,"%s/Z%s_histogram",FN, sufix);
      a.Print_Histogram(fn, Zs, false);
      printf("printed Z histogram!\n");
    
      //------------add to the 3D histogram files--------//
      double* x;
      double* P;

      int Nbins = 12;

      //---taus, both spline and pade
      a.Get_Histogram(Zs, Nbins, x, P, false);

      ZhFile = fopen(ZhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ZhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(ZhFile,"\n");
      fclose(ZhFile);   

      char ZhWFN[300];
      sprintf(ZhWFN, "Z%s_histogram_vs_U.W%.3f.T%.3f",sufix,W,T);
      FILE* ZhWFile = fopen(ZhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ZhWFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fclose(ZhWFile);

      delete [] x;
      delete [] P;

      printf("added to Z_histogram_vs_W and U!\n");

      delete [] mus;
      delete [] Zs;    

    }

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char ZhWFN[300];
      sprintf(ZhWFN, "Z%s_histogram_vs_U.W%.3f.T%.3f",sufix,W,T);
      FILE* ZhWFile = fopen(ZhWFN,"a");
      fprintf(ZhWFile,"\n");
      fclose(ZhWFile);
    }

  }
  }

  return 0;
}
*/


//============================================== tau histogram from real-axis self-consistency ================================================//
/*int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );  

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char tauhWFN[300];
      sprintf(tauhWFN, "tau_RASC_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"w");
      fclose(tauhWFile);
    }
    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { //prepare files for which names need to be prepared only once
      char tauhFN[300];
      sprintf(tauhFN, "tau_RASC_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* tauhFile = fopen(tauhFN,"w");
      fclose(tauhFile);
 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      a.ReadFromUberMinimalFiles(FN);
      double* mus = a.Get_mus();

      complex<double>* Delta;
      double* iw;
      int Miw;
      sprintf(fn,"%s/Delta",FN);
      if (not ReadFunc(fn, Miw, Delta, iw, false)) continue;
      //enforce PH symmetry
      printf("------- read Delta!\n\n");

      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);

      printf(">>>>>>> preparing Deltaw0 and Gw0, pade and spline...\n");

      //spline
      double ImDeltaw0 = imag(CubicFrom4points(Delta,  iw));      

      //pade
      double ImDeltaw0_pade = imag( pade( 200, iw, Delta, 0.0 ) );      
     
      complex<double>* Gw0s = new complex<double>[Nsites];
      complex<double>* Gw0s_pade = new complex<double>[Nsites];

      int Na = a.GetActualN();
      for(int i=0; i<Nsites; i++)
      {
        //spline
        Gw0s[i] = CubicFrom4points(a.r[i].G,  a.r[i].omega);      

        //pade
        Gw0s_pade[i] =  pade( ((Na>200)?200:Na), a.r[i].omega, a.r[i].G, 0.0 ) ;      
      }
      printf(">>>>>>> DONE!\n");

      //--- taus RASC, spline and pade
     
      double* taus_RASC = new double[Nsites];
      double* taus_RASC_pade = new double[Nsites];

      for(int i=0; i<Nsites; i++)
      { taus_RASC[i] = ImDeltaw0 - imag(Gw0s[i])/sqr(abs(Gw0s[i]));
        taus_RASC_pade[i] = ImDeltaw0_pade - imag(Gw0s_pade[i])/sqr(abs(Gw0s_pade[i]));  
      }

      sprintf(fn,"%s/taus_RASC",FN);
      PrintFunc(fn, Nsites, taus_RASC, mus);
      printf("printed taus RASC!\n");
      
      sprintf(fn,"%s/taus_RASC_pade",FN);
      PrintFunc(fn, Nsites, taus_RASC_pade, mus);
      printf("printed taus RASC!\n");

      sprintf(fn,"%s/tau_RASC_histogram",FN);
      a.Print_Histogram(fn, taus_RASC, true);
      printf("printed tau histogram!\n");

      sprintf(fn,"%s/tau_RASC_pade_histogram",FN);
      a.Print_Histogram(fn, taus_RASC_pade, true);
      printf("printed tau histogram!\n");
    
      //------------add to the 3D histogram files--------//
      double* x;
      double* P;
      double* pade_x;
      double* pade_P;

      int Nbins = 12;

      //---taus, both spline and pade
      a.Get_Histogram(taus_RASC, Nbins, x, P, true);
      a.Get_Histogram(taus_RASC_pade, Nbins, pade_x, pade_P, true);

      tauhFile = fopen(tauhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhFile,"%.15le %.15le %.15le %.15le %.15le\n",W,x[b],P[b],pade_x[b],pade_P[b]);
      fprintf(tauhFile,"\n");
      fclose(tauhFile);   

      char tauhWFN[300];
      sprintf(tauhWFN, "tau_RASC_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhWFile,"%.15le %.15le %.15le %.15le %.15le\n",U,x[b],P[b],pade_x[b],pade_P[b]);
      fclose(tauhWFile);

      delete [] x;
      delete [] P;
      delete [] pade_x;
      delete [] pade_P;

      printf("added to tau_histogram_vs_W and U!\n");

      delete [] iw;
      delete [] mus;
      delete [] Delta;
      delete [] Gw0s;
      delete [] Gw0s_pade;
      delete [] taus_RASC;    
      delete [] taus_RASC_pade;    

    }

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char tauhWFN[300];
      sprintf(tauhWFN, "tau_RASC_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"a");
      fprintf(tauhWFile,"\n");
      fclose(tauhWFile);
    }

  }
  }

  return 0;
}
*/

/*
//================================= measure gap and calculate resistivity in the insulator ==================================//

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    // prepare files
    char rhoTFN[300];
    sprintf(rhoTFN,"rho_MI.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);
    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { //prepare files for which names need to be prepared only once
      
 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      
      //---------------- Gavg and Sigma effective----------------------//
      complex<double>* Gw;
      double* w;
     
      int Mw;
      sprintf(fn,"%s/G_avg_w",FN);
      if (not ReadFunc(fn, Mw, Gw, w, false)) continue;
      //enforce PH symmetry
      printf("------- read G_avg_w!\n\n");
            
      if ( abs(imag(Gw[Mw/2]))/pi > 0.01 ) continue; //if metal don't calculate anything

      int i;
      for(i=Mw/2; i<Mw; i++)
        if (abs(imag(Gw[i]))/pi>0.01) break;
      double E_gap = w[i]-w[Mw/2];
             
      rhoTFile = fopen(rhoTFN,"a");
      fprintf(rhoTFile,"%.15le %.15le %.15le %.15le\n", U, W, exp(E_gap/T), E_gap);
      fclose(rhoTFile);    

    }
    rhoTFile = fopen(rhoTFN,"a");
    fprintf(rhoTFile,"\n");
    fclose(rhoTFile);
  }
  }

  return 0;
}
*/

//================================= pade Delta (and G), calc Sigma_eff from RASC and calc resistivity==================================//
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    // prepare files
    char rhoTFN[300];
    sprintf(rhoTFN,"rho_RASC.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);
    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { //prepare files for which names need to be prepared only once
      
 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;
      
      //---------------- Gavg and Sigma effective----------------------//
      complex<double>* Delta;
      complex<double>* G;
      double* iw;
     
      int M;
      sprintf(fn,"%s/Delta",FN);
      if (not ReadFunc(fn, M, Delta, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) Delta[i] = complex<double>(0.0,imag(Delta[i]));
      sprintf(fn,"%s/Delta_read",FN);
      PrintFunc(fn, M, Delta, iw);
      printf("------- read Delta!\n\n");

      sprintf(fn,"%s/G_avg",FN);
      if (not ReadFunc(fn, M, G, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) G[i] = complex<double>(0.0,imag(G[i]));
      printf("------- read G_avg!\n\n");
      
      int Mw = 2000;
      complex<double>* Gw = new complex<double>[Mw];
      complex<double>* Sigmaw = new complex<double>[Mw];
      complex<double>* Deltaw = new complex<double>[Mw];
      
      double* w = new double[Mw];
      double wmax = max(U/2.0,W/2.0) + 2.0;
      for(int i=0; i<Mw; i++)
        w[i] = - wmax + i * 2.0*wmax/(Mw-1.0);     

      //PrintFunc("w_created",Mw,w); 
      printf("------- about to do pade...\n\n");

      pade( 2000, iw, Delta, 
            Mw, w, Deltaw );
      sprintf(fn,"%s/Delta_w",FN);
      PrintFunc(fn, Mw, Deltaw, w);
      
      //pade G_avg ...
      pade( 2000, iw, G, 
            Mw, w, Gw );
      sprintf(fn,"%s/G_avg_w",FN);
      PrintFunc(fn, Mw, Gw, w);

      // ... or read it if you have it
      //  sprintf(fn,"%s/G_avg_w",FN);
      //  if (not ReadFunc(fn, Mw, Gw, w, false)) continue;

      printf(">>>>>>>> Did pade on Delta and G and printed to files!\n\n");           


      complex<double>* Sigma_eff_w_RASC = new complex<double>[Mw];
      for(int i=0; i<Mw; i++)
        Sigma_eff_w_RASC[i] = w[i]+U/2.0-Deltaw[i]-1.0/Gw[i];

      sprintf(fn,"%s/Sigma_eff_w_RASC",FN);
      PrintFunc(fn, Mw, Sigma_eff_w_RASC, w);

      printf(">>>>>>>> Printed Sigma_eff_w obtained from self-consistency on the real axis!\n\n");

      //-------- conductivity from real axis Sigma effective -----------//      

      printf("------- resistivity...\n");
      double sgm =  CubicConductivity(Mw, Sigma_eff_w_RASC, w, U/2.0, T, 16, 400, 1000.0);
      double sgm2 =  CubicConductivity(Mw, Sigma_eff_w_RASC, w, U/2.0, T, 16, 400, 1000.0, true);

      rhoTFile = fopen(rhoTFN,"a");
      fprintf(rhoTFile,"%.15le %.15le %.15le %.15le\n", U, W, 1.0/sgm, 1.0/sgm2);
      fclose(rhoTFile);

      printf("------- resistivity done!\n");

      delete [] Delta;
      delete [] Deltaw;
      delete [] Gw;
      delete [] Sigma_eff_w_RASC;
      delete [] iw;
      delete [] w;

    }
    rhoTFile = fopen(rhoTFN,"a");
    fprintf(rhoTFile,"\n");
    fclose(rhoTFile);
  }
  }

  return 0;
}

*/

//=================================== ALL ============================================//
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    // prepare files
    
    char rhoTFN[300];
    sprintf(rhoTFN,"rho.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);
    

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char Aw0sWFN[300];
      sprintf(Aw0sWFN, "Aw0s_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0sWFile = fopen(Aw0sWFN,"w");
      fclose(Aw0sWFile);

      char Aw0hWFN[300];
      sprintf(Aw0hWFN, "Aw0_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0hWFile = fopen(Aw0hWFN,"w");
      fclose(Aw0hWFile);

      char nhWFN[300];
      sprintf(nhWFN, "n_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* nhWFile = fopen(nhWFN,"w");
      fclose(nhWFile); 

      char tauhWFN[300];
      sprintf(tauhWFN, "tau_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"w");
      fclose(tauhWFile);

      char rmuhWFN[300];
      sprintf(rmuhWFN, "rmu_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* rmuhWFile = fopen(rmuhWFN,"w");
      fclose(rmuhWFile); 

      char ReSigw0hWFN[300];
      sprintf(ReSigw0hWFN, "ReSigw0_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ReSigw0hWFile = fopen(ReSigw0hWFN,"w");
      fclose(ReSigw0hWFile); 
    }
    
    // ===================== ITERATE ======================== //
    for(double U = 0.25; U<9.6; U+=0.25)
    { //prepare files for which names need to be prepared only once
      char Aw0hFN[300];
      sprintf(Aw0hFN, "Aw0_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* Aw0hFile = fopen(Aw0hFN,"w");
      fclose(Aw0hFile);

      char Aw0sUFN[300];
      sprintf(Aw0sUFN, "Aw0s_vs_W.U%.3f.T%.3f",U,T);
      FILE* Aw0sUFile = fopen(Aw0sUFN,"w");
      fclose(Aw0sUFile);

      char nhFN[300];
      sprintf(nhFN, "n_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* nhFile = fopen(nhFN,"w");
      fclose(nhFile); 

      char tauhFN[300];
      sprintf(tauhFN, "tau_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* tauhFile = fopen(tauhFN,"w");
      fclose(tauhFile);

      char rmuhFN[300];
      sprintf(rmuhFN, "rmu_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* rmuhFile = fopen(rmuhFN,"w");
      fclose(rmuhFile);

      char ReSigw0hFN[300];
      sprintf(ReSigw0hFN, "ReSigw0_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* ReSigw0hFile = fopen(ReSigw0hFN,"w");
      fclose(ReSigw0hFile);
 
    for(double W = 0.25; W<12.1; W+=0.25)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      char cmd[300];
      sprintf(cmd,"rm %s/mu*",FN);
      system(cmd);
      
      a.ReadFromUberMinimalFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);
     
      //prints ns and ns_sorted
      a.Print_ns_and_mus(FN);
      printf("printed ns and mus!\n");

      
      //--- SPLINE

      double* Aw0s = a.Get_Aw0s();
      double* taus = a.Get_taus();  
      double* mus = a.Get_mus();
      double* rmus = a.Get_renormalized_mus();
      double* ReSigw0s = a.Get_ReSigw0s();

      //print sorted
      sprintf(fn,"%s/Aw0s",FN);
      PrintFunc(fn, Nsites, Aw0s, mus);

      Aw0sUFile = fopen(Aw0sUFN,"a");
      for(int i=0; i<Nsites; i++) 
        fprintf(Aw0sUFile,"%.15le %.15le %.15le\n", W, mus[i], Aw0s[i]);
      fprintf(Aw0sUFile,"\n"); 
      fclose(Aw0sUFile);
      
      char Aw0sWFN[300];
      sprintf(Aw0sWFN, "Aw0s_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0sWFile = fopen(Aw0sWFN,"a");
      for(int i=0; i<Nsites; i++) 
        fprintf(Aw0sWFile,"%.15le %.15le %.15le\n", U, mus[i], Aw0s[i]);
      fclose(Aw0sWFile);         
      printf("printed Aw0s!\n");

      sprintf(fn,"%s/taus",FN);
      PrintFunc(fn, Nsites, taus, mus);
      printf("printed taus!\n");
      
      //print histograms
      sprintf(fn,"%s/Aw0_histogram",FN);
      a.Print_Histogram(fn, Aw0s);
      printf("printed Aw0 histogram!\n");

      sprintf(fn,"%s/tau_histogram",FN);
      a.Print_Histogram(fn, taus, true);
      printf("printed tau histogram!\n");

      sprintf(fn,"%s/rmu_histogram",FN);
      a.Print_Histogram(fn, rmus);
      printf("printed renormalized mu histogram!\n");

      sprintf(fn,"%s/ReSigw0_histogram",FN);
      a.Print_Histogram(fn, ReSigw0s, true);
      printf("printed ReSigw0 histogram!\n");


      //--- PADE

      double* pade_Aw0s = a.Get_pade_Aw0s();
      double* pade_taus = a.Get_pade_taus();  

      //print sorted
      sprintf(fn,"%s/Aw0s_pade",FN);
      PrintFunc(fn, Nsites, pade_Aw0s, mus);
      printf("printed Aw0s!\n");

      sprintf(fn,"%s/taus_pade",FN);
      PrintFunc(fn, Nsites, pade_taus, mus);
      printf("printed taus!\n");
      
      //print histograms
      sprintf(fn,"%s/Aw0_pade_histogram",FN);
      a.Print_Histogram(fn, pade_Aw0s);
      printf("printed Aw0 histogram!\n");

      sprintf(fn,"%s/tau_pade_histogram",FN);
      a.Print_Histogram(fn, pade_taus, true);
      printf("printed tau histogram!\n");

      double* ns = a.Get_ns();     
      sprintf(fn,"%s/n_histogram",FN);
      a.Print_Histogram(fn, ns);
      printf("printed n histogram!\n");
     
      //------------add to the 3D histogram files--------//
      double* x;
      double* P;
      double* pade_x;
      double* pade_P;

      int Nbins = 12;

      //---Aw0s, both spline and pade
      a.Get_Histogram(Aw0s, Nbins, x, P);
      a.Get_Histogram(pade_Aw0s, Nbins, pade_x, pade_P);

      Aw0hFile = fopen(Aw0hFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(Aw0hFile,"%.15le %.15le %.15le %.15le %.15le\n",W,x[b],P[b],pade_x[b],pade_P[b]);
      fprintf(Aw0hFile,"\n");
      fclose(Aw0hFile);   

      char Aw0hWFN[300];
      sprintf(Aw0hWFN, "Aw0_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0hWFile = fopen(Aw0hWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(Aw0hWFile,"%.15le %.15le %.15le %.15le %.15le\n",U,x[b],P[b],pade_x[b],pade_P[b]);
      fclose(Aw0hWFile);  

      delete [] x;
      delete [] P;
      delete [] pade_x;
      delete [] pade_P;
      printf("added to Aw0_histogram_vs_W and U!\n");

      //---taus, both spline and pade
      a.Get_Histogram(taus, Nbins, x, P, true);
      a.Get_Histogram(pade_taus, Nbins, pade_x, pade_P, true);

      tauhFile = fopen(tauhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhFile,"%.15le %.15le %.15le %.15le %.15le\n",W,x[b],P[b],pade_x[b],pade_P[b]);
      fprintf(tauhFile,"\n");
      fclose(tauhFile);   

      char tauhWFN[300];
      sprintf(tauhWFN, "tau_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhWFile,"%.15le %.15le %.15le %.15le %.15le\n",U,x[b],P[b],pade_x[b],pade_P[b]);
      fclose(tauhWFile);

      delete [] x;
      delete [] P;
      delete [] pade_x;
      delete [] pade_P;

      printf("added to tau_histogram_vs_W!\n");

      //---ns
      a.Get_Histogram(ns, Nbins, x, P);

      nhFile = fopen(nhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(nhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(nhFile,"\n");
      fclose(nhFile);   

      char nhWFN[300];
      sprintf(nhWFN, "n_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* nhWFile = fopen(nhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(nhWFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fclose(nhWFile); 

      delete [] x;
      delete [] P;
      printf("added to n_histogram_vs_W!\n");

      //---Renomralized mus
      a.Get_Histogram(rmus, Nbins, x, P);

      rmuhFile = fopen(rmuhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(rmuhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(rmuhFile,"\n");
      fclose(rmuhFile);   

      char rmuhWFN[300];
      sprintf(rmuhWFN, "rmu_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* rmuhWFile = fopen(rmuhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(rmuhWFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fclose(rmuhWFile); 

      delete [] x;
      delete [] P;
      printf("added to rmu_histogram_vs_W!\n");

      //---ReSigw0s
      a.Get_Histogram(ReSigw0s, Nbins, x, P);

      ReSigw0hFile = fopen(ReSigw0hFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ReSigw0hFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(ReSigw0hFile,"\n");
      fclose(ReSigw0hFile);   

      char ReSigw0hWFN[300];
      sprintf(ReSigw0hWFN, "ReSigw_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ReSigw0hWFile = fopen(ReSigw0hWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ReSigw0hWFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fclose(ReSigw0hWFile); 

      delete [] x;
      delete [] P;
      printf("added to rmu_histogram_vs_W!\n");

      delete [] ns;
      delete [] mus;
      delete [] Aw0s;
      delete [] taus;    
      delete [] pade_Aw0s;
      delete [] pade_taus;    

      //---------------- Gavg and Sigma effective----------------------//
      complex<double>* Sigma;
      complex<double>* G;
      double* iw;
     
      int M;
      sprintf(fn,"%s/Sigma_eff",FN);
      if (not ReadFunc(fn, M, Sigma, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) Sigma[i] = complex<double>(U/2.0,imag(Sigma[i]));
      //sprintf(fn,"%s/Sigma_eff_read",FN);
      //PrintFunc(fn, M, Sigma, iw);
      printf("------- read Sigma effective!\n\n");
      
      sprintf(fn,"%s/G_avg",FN);
      if (not ReadFunc(fn, M, G, iw, false)) continue;
      //enforce PH symmetry
      for(int i=0; i<M; i++) G[i] = complex<double>(0.0,imag(G[i]));
      //sprintf(fn,"%s/G_avg_read",FN);
      //PrintFunc(fn, M, G, iw);
      printf("------- read G average!\n\n");

      int Mw = 2000;
      complex<double>* Gw = new complex<double>[Mw];
      complex<double>* Sigmaw = new complex<double>[Mw];
      double* w = new double[Mw];
      double wmax = max(U/2.0,W/2.0) + 2.0;
      for(int i=0; i<Mw; i++)
        w[i] = - wmax + i * 2.0*wmax/(Mw-1.0);     

      //PrintFunc("w_created",Mw,w); 
      printf("------- about to do pade...\n\n");

     //this is already done, we don't need it
   //   pade( 2000, iw, G, 
   //         Mw, w, Gw );
   //   sprintf(fn,"%s/G_avg_w",FN);
   //   PrintFunc(fn, Mw, Gw, w);
   
      pade( 2000, iw, Sigma, 
            Mw, w, Sigmaw );
      sprintf(fn,"%s/Sigma_eff_w",FN);
      PrintFunc(fn, Mw, Sigmaw, w);

      printf(">>>>>>>> Did pade on G_avg and Sigma_eff and printed to files!\n\n");

      //-------- conductivity from imag axis Sigma effective -----------//      
      printf("------- sigma(nu)...\n");
      complex<double> Lambda0 = Lambda(M, Sigma, iw, U/2.0, 16, 0); //calculate zero frequency current-current correlation only once

      int Nn = 400;
      complex<double>* sigma = new complex<double>[Nn];
      double* nu = new double[Nn]; 
      for(int n=0; n<Nn; n++)
      {  sigma[n] = OpticalConductivity(M, Sigma, iw, U/2.0, 16, n+1, Lambda0);
         nu[n] = 2.0*(n+1)*pi*T;
      }
      
      sprintf(fn,"%s/sigma_nu",FN);
      PrintFunc(fn,Nn,sigma,nu);

      double rho = 1.0/real(pade( Nn, nu, sigma, 0.0 ));
      printf("------- sigma(nu) DONE!\n");

      //-------- conductivity from real axis Sigma effective -----------//      

      printf("------- resistivity...\n");
      double sgm =  CubicConductivity(Mw, Sigmaw, w, U/2.0, T, 16, 400, 1000.0);

      rhoTFile = fopen(rhoTFN,"a");
      fprintf(rhoTFile,"%.15le %.15le %.15le %.15le\n", U, W, 1.0/sgm, rho);
      fclose(rhoTFile);

      printf("------- resistivity done!\n");

      delete [] G;
      delete [] Sigma;
      delete [] Gw;
      delete [] Sigmaw;
      delete [] iw;
      delete [] w;
      delete [] nu;

      //release H0     
    }
    rhoTFile = fopen(rhoTFN,"a");
    fprintf(rhoTFile,"\n");
    fclose(rhoTFile);

    for(double W = 0.25; W<12.1; W+=0.25)
    {
      char Aw0sWFN[300];
      sprintf(Aw0sWFN, "Aw0s_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0sWFile = fopen(Aw0sWFN,"a");
      fprintf(Aw0sWFile,"\n");
      fclose(Aw0sWFile);

      char Aw0hWFN[300];
      sprintf(Aw0hWFN, "Aw0_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* Aw0hWFile = fopen(Aw0hWFN,"a");
      fprintf(Aw0hWFile,"\n");
      fclose(Aw0hWFile);  

      char tauhWFN[300];
      sprintf(tauhWFN, "tau_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhWFile = fopen(tauhWFN,"a");
      fprintf(tauhWFile,"\n");
      fclose(tauhWFile);

      char nhWFN[300];
      sprintf(nhWFN, "n_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* nhWFile = fopen(nhWFN,"a");
      fprintf(nhWFile,"\n");
      fclose(nhWFile); 

      char rmuhWFN[300];
      sprintf(rmuhWFN, "rmu_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* rmuhWFile = fopen(rmuhWFN,"a");
      fprintf(rmuhWFile,"\n");
      fclose(rmuhWFile); 

      char ReSigw0hWFN[300];
      sprintf(ReSigw0hWFN, "ReSigw0_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ReSigw0hWFile = fopen(ReSigw0hWFN,"a");
      fprintf(ReSigw0hWFile,"\n");
      fclose(ReSigw0hWFile); 
    }

  }
  }

  return 0;
}
*/

//================================== fixedW run =============================================//
/*int main(int argc, char* argv [])
{
//  if(argc<4) exit(0);
//  double W = atof(argv[1]);
//  double T = atof(argv[2]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 100; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  for(double T=0.01; T<0.011; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );
    
    for(double W = 2.0; W<7.1; W+=1.0)
    { //if ((W>4.8)and(W<5.2)) continue;
      char nhFN[300];
      sprintf(nhFN, "n_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* nhFile = fopen(nhFN,"w");
      fclose(nhFile); 

      char tauhFN[300];
      sprintf(tauhFN, "tau_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhFile = fopen(tauhFN,"w");
      fclose(tauhFile);
 
    for(double U = 0.6; U<W+2.1; U+=0.2)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IACPA.cubic.averageG.W%.3f.T%.3f.U%.3f",W,T,U);
      
      char cmd[300];
      sprintf(cmd,"rm %s/mu*",FN);
      system(cmd);

      a.ReadFromFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);
     
      //prints ns and ns_sorted
      a.Print_ns_and_mus(FN);
      printf("printed ns and mus!\n");

      //print sorted
      char fn[300];
      double* Aw0s = a.Get_Aw0s();
      double* taus = a.Get_taus();  
      double* mus = a.Get_mus();
  
      sprintf(fn,"%s/Aw0s",FN);
      PrintFunc(fn, Nsites, Aw0s, mus);
      printf("printed Aw0s sorted!\n");

      sprintf(fn,"%s/taus",FN);
      PrintFunc(fn, Nsites, taus, mus);
      printf("printed sorted taus!\n");
      
      //print histograms
      sprintf(fn,"%s/tau_histogram",FN);
      a.Print_Histogram(fn, taus, true);
      printf("printed tau histogram!\n");

      double* ns = a.GetSorted_ns();     
      sprintf(fn,"%s/n_histogram",FN);
      a.Print_Histogram(fn, ns);
      printf("printed n histogram!\n");
     
      //add to the 3D histogram files
      double* x;
      double* P;

      int Nbins = 12;

      a.Get_Histogram(ns, Nbins, x, P);

      nhFile = fopen(nhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(nhFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fprintf(nhFile,"\n");
      fclose(nhFile);   
      delete [] x;
      delete [] P;
      printf("added to n_histogram_vs_U!\n");

      a.Get_Histogram(taus, Nbins, x, P, true);

      tauhFile = fopen(tauhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(tauhFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fprintf(tauhFile,"\n");
      fclose(tauhFile);   
      delete [] x;
      delete [] P;
      printf("added to tau_histogram_vs_U!\n");

      delete [] ns;
      delete [] mus;
      delete [] Aw0s;
      delete [] taus;    

      //Pade on Green's function and Self-energy
      //for(int id=0; id<Nsites; id++)
      //{ char pFN[300];
      //  sprintf(pFN,"%sGw.%d",FN,id);
      //  PadeToFile( 2000, a.r[id].G,  a.r[id].omega, pFN, 500, 3.0 );
      //  sprintf(pFN,"%sSigw.%d",FN,id);
      //  PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 500, 3.0 );

      //  printf("pade, site id=%d : DONE!\n",id);
      //}
      
      //release H0     
    }
  }
  }

  return 0;
}*/
