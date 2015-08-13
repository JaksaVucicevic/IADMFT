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


//============================================== Pade G and Delta, get Sigma from RASC, calc local resistivity and histogram================================================//


int main(int argc, char* argv [])
{
  if(argc<4) exit(0);
  double W = atof(argv[1]);
  double T = atof(argv[2]);
  double U = atof(argv[3]);
//  double mu=U/2.0; 

   // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  //double t=1.0/6.0; //moved to global
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];
  char pFN[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.002; T<0.52; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );
    
    //for(double W=5.0; W<7.1; W+=100.0)
    { 
      
    //for(double U=W-1.0; U<W+1.0; U+=0.05)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      
      bool successful = a.ReadFromMinimalFiles(FN);
      if (not successful) return -1;
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);    
       
      int actualN = a.GetActualN();

      double* rhos = new double[Nsites];
      for(int id=0; id<Nsites; id++)
      { 
        double wmax = 7.0;
        int NwG = 1000;
        double* wG = new double[NwG]; //omega grid
        complex<double>* Gw = new complex<double>[NwG]; //green's function
        complex<double>* Deltaw = new complex<double>[NwG]; //hybridization function
        for(int i=0; i<NwG; i++)
          wG[i] = -wmax + i*2.0*wmax/((double)NwG-1.0);
          
        pade( actualN, a.r[id].omega, a.r[id].G, 
             NwG, wG, Gw );
        pade( actualN, a.r[id].omega, a.r[id].Delta, 
             NwG, wG, Deltaw );

        sprintf(pFN,"%s/Gw.%d",FN,id);
        PrintFunc(pFN, NwG, Gw, wG); 
        sprintf(pFN,"%s/Deltaw.%d",FN,id);
        PrintFunc(pFN, NwG, Deltaw, wG);
        
        int NwSigma = 10000;
        double* wSigma = new double[NwSigma];
        complex<double>* Sigmaw = new complex<double>[NwSigma];
        
        for(int i=0; i<NwSigma; i++)
        { wSigma[i] = -wmax + i*2.0*wmax/((double)NwSigma-1.0);
          Sigmaw[i] = wSigma[i] + a.r[id].mu - interpl(NwG, Deltaw, wG, wSigma[i]) - 1.0/interpl(NwG, Gw, wG, wSigma[i]);
        }

        sprintf(pFN,"%s/Sigmaw.%d",FN,id);
        PrintFunc(pFN,NwSigma,Sigmaw,wSigma);

        rhos[id] = 1.0/CubicConductivity(NwSigma, Sigmaw, wSigma, a.r[id].mu, T, 16, 400, 1000.0, false);

        delete [] wG;
        delete [] Gw;
        delete [] Deltaw;        
        delete [] wSigma;
        delete [] Sigmaw;

        printf("pade, site id=%d : DONE!\n",id);
      }
      double* mus = a.Get_mus();
      double* mus_sorted = a.GetSorted(mus);
      double* rhos_sorted = a.GetSorted(rhos);
      double* ns = a.Get_ns();
      double* ns_sorted = a.GetSorted(ns);
      double* Aw0s = a.Get_Aw0s();
      double* Aw0s_pade = a.Get_pade_Aw0s();
      double* Aw0s_sorted = a.GetSorted(Aw0s);
      double* Aw0s_pade_sorted = a.GetSorted(Aw0s_pade);
      double* Zs_RASC = a.Get_Zs_RASC();
      double* taus_RASC = a.Get_taus_RASC();
      double* Zs_RASC_sorted = a.GetSorted(Zs_RASC);
      double* taus_RASC_sorted = a.GetSorted(taus_RASC);

     
      int Nbins = 15;
      double* x = new double[Nbins];
      double* P = new double[Nbins];
      Histogram(Nsites, rhos, Nbins, x, P, false);

      double* logx = new double[Nbins];
      double* logP = new double[Nbins];
      Histogram(Nsites, rhos, Nbins, logx, logP, true);

      sprintf(pFN,"%s/Aw0s",FN);
      PrintFunc(pFN, Nsites, Aw0s, mus);

      sprintf(pFN,"%s/Aw0s_sorted",FN);
      PrintFunc(pFN, Nsites, Aw0s_sorted, mus_sorted);

      sprintf(pFN,"%s/Aw0s_pade",FN);
      PrintFunc(pFN, Nsites, Aw0s_pade, mus);

      sprintf(pFN,"%s/Aw0s_pade_sorted",FN);
      PrintFunc(pFN, Nsites, Aw0s_pade_sorted, mus_sorted);

      sprintf(pFN,"%s/ns",FN);
      PrintFunc(pFN, Nsites, ns, mus);

      sprintf(pFN,"%s/ns_sorted",FN);
      PrintFunc(pFN, Nsites, ns_sorted, mus_sorted);

      sprintf(pFN,"%s/Zs_RASC",FN);
      PrintFunc(pFN, Nsites, Zs_RASC, mus);

      sprintf(pFN,"%s/Zs_RASC_sorted",FN);
      PrintFunc(pFN, Nsites, Zs_RASC_sorted, mus_sorted);

      sprintf(pFN,"%s/taus_RASC",FN);
      PrintFunc(pFN, Nsites, taus_RASC, mus);

      sprintf(pFN,"%s/taus_RASC_sorted",FN);
      PrintFunc(pFN, Nsites, taus_RASC_sorted, mus_sorted);

      sprintf(pFN,"%s/rhos",FN);
      PrintFunc(pFN, Nsites, rhos, mus);

      sprintf(pFN,"%s/rhos_sorted",FN);
      PrintFunc(pFN, Nsites, rhos_sorted, mus_sorted);
      
      sprintf(pFN,"%s/rho_histogram",FN);
      PrintFunc(pFN, Nbins, P, x);

      sprintf(pFN,"%s/rho_log_histogram",FN);
      PrintFunc(pFN, Nbins, logP, logx);


      delete [] rhos;
      delete [] rhos_sorted;

      delete [] P;
      delete [] logP;
      delete [] x;
      delete [] logx;

      delete [] mus;
      delete [] mus_sorted;
      delete [] rhos_sorted;
      delete [] ns;
      delete [] ns_sorted;
      delete [] Aw0s;
      delete [] Aw0s_pade;
      delete [] Aw0s_sorted;
      delete [] Aw0s_pade_sorted;
      delete [] Zs_RASC;
      delete [] taus_RASC;
      delete [] Zs_RASC_sorted;
      delete [] taus_RASC_sorted;
    }
  }
  }

  return 0;
}




/*

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 

  char sufix[300];
  //sprintf(sufix,"_RASC");
  sprintf(sufix,"");

  int Nsites = 216; //Nepsilons
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  //for(double T=0.1; T>0.0009; T/=10.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );  

    for(double W = 2.0; W<7.0; W+=1.0)
    {
      char ZhWFN[300];
      sprintf(ZhWFN, "Z%s_histogram_vs_U.W%.3f.T%.3f",sufix,W,T);
      FILE* ZhWFile = fopen(ZhWFN,"w");
      fclose(ZhWFile);
    }
    
    // ===================== ITERATE ======================== //
    for(double U = 0.6; U<9.0; U+=0.1)
    { //prepare files for which names need to be prepared only once
      char ZhFN[300];
      sprintf(ZhFN, "Z%s_histogram_vs_W.U%.3f.T%.3f",sufix,U,T);
      FILE* ZhFile = fopen(ZhFN,"w");
      fclose(ZhFile);
 
    for(double W = ceil(max(1.999,U-2.0)); W<7.0; W+=1.0)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      a.ReadFromFiles(FN);
      double* mus = a.Get_mus();

      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);

      double* Zs = a.Get_Zs_RASC();
      a.GetSorted(Zs);

      sprintf(fn,"%s/Zs%s",FN,sufix);
      PrintFunc(fn, Nsites, Zs, mus);
      printf("printed Zs!\n");
      
      sprintf(fn,"%s/Z%s_histogram",FN,sufix);
      a.Print_Histogram(fn, Zs, false);
      printf("printed Z%s histogram!\n",sufix);
    
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

    for(double W = 2.0; W<7.0; W+=1.0)
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
/*
int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
  double T = atof(argv[1]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  int Nsites = 216; //Nsites
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

    for(double W = 2.0; W<7.0; W+=1.0)
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

      char ZhWFN[300];
      sprintf(ZhWFN, "Z_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ZhWFile = fopen(ZhWFN,"w");
      fclose(ZhWFile);
    }
    
    // ===================== ITERATE ======================== //
    for(double U = 0.6; U<9.0; U+=0.1)
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

      char ZhFN[300];
      sprintf(ZhFN, "Z_histogram_vs_W.U%.3f.T%.3f",U,T);
      FILE* ZhFile = fopen(ZhFN,"w");
      fclose(ZhFile);
 
    for(double W = ceil(max(1.999,U-2.0)); W<7.0; W+=1.0)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      char fn[300];
      sprintf(fn,"%s/0",FN);
      if (not FileExists(fn)) continue;

      //char cmd[300];
      //sprintf(cmd,"rm %s/mu*",FN);
      //system(cmd);
      
      a.ReadFromFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);
     
      //prints ns and ns_sorted
      a.Print_ns_and_mus(FN);
      printf("printed ns and mus!\n");

      
      //--- SPLINE

      double* Aw0s = a.Get_Aw0s();
      a.GetSorted(Aw0s);
      double* taus = a.Get_taus();  
      a.GetSorted(taus);
      double* mus = a.Get_mus();
      a.GetSorted(mus);
      double* rmus = a.Get_renormalized_mus();
      a.GetSorted(rmus);
      double* ReSigw0s = a.Get_ReSigw0s();
      a.GetSorted(ReSigw0s);

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

      double* Zs = a.Get_Zs();
      a.GetSorted(Zs);

      sprintf(fn,"%s/Zs",FN);
      PrintFunc(fn, Nsites, Zs, mus);
      printf("printed Zs!\n");
      
      sprintf(fn,"%s/Z_histogram",FN);
      a.Print_Histogram(fn, Zs, false);
      printf("printed Z histogram!\n");


      //--- PADE

      double* pade_Aw0s = a.Get_pade_Aw0s();
      a.GetSorted(pade_Aw0s);
      double* pade_taus = a.Get_pade_taus();  
      a.GetSorted(pade_taus);

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
      a.GetSorted(ns);    
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

      //---Zs
      a.Get_Histogram(Zs, Nbins, x, P, false);

      ZhFile = fopen(ZhFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ZhFile,"%.15le %.15le %.15le\n",W,x[b],P[b]);
      fprintf(ZhFile,"\n");
      fclose(ZhFile);   

      char ZhWFN[300];
      sprintf(ZhWFN, "Z_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ZhWFile = fopen(ZhWFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(ZhWFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fclose(ZhWFile);

      delete [] x;
      delete [] P;
      printf("added to rmu_histogram_vs_W!\n");

      delete [] ns;
      delete [] mus;
      delete [] Aw0s;
      delete [] taus;  
      delete [] Zs;  
      delete [] pade_Aw0s;
      delete [] pade_taus;    
      delete [] rmus;
      delete [] ReSigw0s;

      //release H0     
    }

    for(double W = 1.0; W<7.0; W+=1.0)
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

      char ZhWFN[300];
      sprintf(ZhWFN, "Z_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* ZhWFile = fopen(ZhWFN,"a");
      fprintf(ZhWFile,"\n");
      fclose(ZhWFile);
    }

  }
  }

  return 0;
}
*/

/*
int main(int argc, char* argv [])
{
//  if(argc<4) exit(0);
//  double W = atof(argv[1]);
//  double T = atof(argv[2]);
//  double U = atof(argv[3]);
  double mu=U/2.0; 

   // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=1.0/6.0; 
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  for(double T=0.01; T<0.011; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );
    
    for(double W=5.0; W<7.1; W+=100.0)
    { 
      char Aw0FN[300];
      sprintf(Aw0FN,"Aw0_from_G_avg_pade.W%.3f.T%.3f",W,T);
      FILE* Aw0File = fopen(Aw0FN,"w");
      fclose(Aw0File);

    for(double U=1.0; U<W+2.0; U+=0.1)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      
      a.ReadFromFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);    

      //Pade on Green's function and Self-energy
      complex<double>* G = new complex<double>[Nw];
      complex<double>* Sigma = new complex<double>[Nw];

      a.GetAverageGandSigma(G,Sigma);
  
      complex<double> Gw0 = pade( 2000,  a.r[0].omega, G, 0.0 );
      Aw0File = fopen(Aw0FN,"a"); 
      fprintf(Aw0File,"%.15le %.15le %.15le %.15le\n", U, -(1.0/pi)*imag(Gw0), real(Gw0), imag(Gw0) );
      fclose(Aw0File);
      
      //char pFN[300];
      //sprintf(pFN,"%s/G_avg_w",FN);
      //PadeToFile( 2000, G,  a.r[0].omega, pFN, 1000, 5.0 );
        
      //for(int id=0; id<Nsites; id++)
      //{ 
      //  sprintf(pFN,"%s/Gw.%d",FN,id);
      //  PadeToFile( 2000, a.r[id].G,  a.r[id].omega, pFN, 1000, 5.0 );
      //  //sprintf(pFN,"%s/Sigw.%d",FN,id);
      //  //PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 1000, 5.0 );

      //  printf("pade, site id=%d : DONE!\n",id);
      //}
      
      delete [] G;
      delete [] Sigma;
      
      //release H0     
    }
  }
  }

  return 0;
}


*/















//======================================== all data ================================================//
/*int main(int argc, char* argv [])
{
   // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=1.0/6.0; 
  //--------------- StatDMFT -------------------//

  char FN[300];
  char cmd[300];

  //matsubara grid
  int Nw = pow(2,13);

  for(double T=0.01; T<0.011; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );
    
    for(double W=5.0; W<7.1; W+=100.0)
    { //if ((W>4.8)and(W<5.2)) continue;
      char nhFN[300];
      sprintf(nhFN, "n_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* nhFile = fopen(nhFN,"w");
      fclose(nhFile); 

      char tauhFN[300];
      sprintf(tauhFN, "tau_histogram_vs_U.W%.3f.T%.3f",W,T);
      FILE* tauhFile = fopen(tauhFN,"w");
      fclose(tauhFile);
 
    for(double U=1.0; U<W+2.01; U+=0.1)
    { //==========//

      //load the result and clear previous sorted files
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      
      char cmd[300];
      sprintf(cmd,"rm %s/mu*",FN);
      system(cmd);

      a.ReadFromFiles(FN);
      printf("======================= FILE LOADED: U: %.3f T: %.3f W: %.3f =======================\n",U,T,W);

      //print sorted results  
      a.PrintSortedMinimal(FN, 10.0);
      printf("printed sorted results!\n");     

      //prints ns and ns_sorted
      a.Print_ns_and_mus(FN);
      printf("printed ns and mus!\n");

      //print sorted
      char fn[300];
      double* Aw0s = a.GetSorted_Aw0s();
      double* taus = a.GetSorted_taus();  
      double* mus = a.GetSorted_mus();
  
      sprintf(fn,"%s/Aw0s_sorted",FN);
      PrintFunc(fn, Nsites, Aw0s, mus);
      printf("printed Aw0s sorted!\n");

      sprintf(fn,"%s/taus_sorted",FN);
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

      a.Get_Histogram(taus, Nbins, x, P, true); //!! should be redone on logscale

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

      //print average G and Sigma
      complex<double>* G = new complex<double>[Nw];
      complex<double>* Sigma = new complex<double>[Nw];

      a.GetAverageGandSigma(G,Sigma);
  
      sprintf(fn,"%s/G_avg",FN);
      PrintFunc(fn, Nw, G, a.r[0].omega);
      sprintf(fn,"%s/Sigma_avg",FN);
      PrintFunc(fn, Nw, Sigma, a.r[0].omega);
      printf("Printed average G and Sigma\n");  

      delete [] G;
      delete [] Sigma;
  
    }
  }
  }

  return 0;
}*/
