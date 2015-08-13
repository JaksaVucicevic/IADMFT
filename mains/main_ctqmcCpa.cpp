#include "../source/IAresArray.h"
#include "../source/IAResult.h"
#include "../source/IAGRID.h"
#include "../source/arrayInitializers.h"
#include "../source/ctqmcCPA.h"
#include "../source/ctqmcSIAM.h"
#include "../source/routines.h"
#include "../source/pade.h"
#include <omp.h>

#include <cstdlib>
#include <cstdio>

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double W = atof(argv[1]);
  //double T = atof(argv[1]);
  //double U = atof(argv[2]);

for(double T=0.002; T<0.52; T*=2.0)
{  int Nepsilons = 100; 
  double t=0.5; //for the cubic this makes no difference, except for the initial guess which may be semi-circle

  int Nw = 10000;
  double* NIDOS = new double[Nw];
  double* w = new double[Nw];
  ReadFunc("cubic_dos", Nw, NIDOS, w);

  //reduce the number of points for the integral to work more quickly
  int ratio = 10;
  Nw /= ratio;
  for(int i=0; i<Nw; i++)
  { NIDOS[i] = NIDOS[ratio*i];
    w[i] = w[ratio*i]; 
  }
  
  PrintFunc("cubic_dos_read", Nw, NIDOS, w);

  char FN[300];
  char cmd[300];

  //matsubara grid
  int N = pow(2,13);

  IAGRID g( N, N, T );

  //all the storage arrays
  IAresArray a( Nepsilons, &g );

  //initialized Delta to metallic
  /*for(int i=0; i<Nepsilons; i++)
  for(int j=0; j<N; j++)
    a.r[i].Delta[j] = SemiCircleG(ii*a.r[i].omega[j],t);
  */
  //a.r[0].PrintResult("proba");

  //-------------//
  ctqmcCPA cpa;
  cpa.UseBethe = false; //if true, UseAveragedSigma makes no difference
  cpa.UseAverageSigma = false;
  cpa.NIDOS = NIDOS;
  cpa.w = w;
  cpa.Nw = Nw;
  cpa.PHSymmetry = true;
  cpa.PatchDelta = true;
  cpa.AtomicCutoff = 30.0;

  cpa.SetBroydenOptions(false, false, 0);

  //init ffts
  cpa.fft = new FFT[Nepsilons];
  for(int id=0; id<Nepsilons; id++)
    cpa.fft[id].Initialize(N, T, a.r[id].omega, a.r[id].tau);   


char prefix[300] = "IACPA.cubic.averageG";

//for(double W = 0.25; W<12.1; W+=0.25)
{    
for(double U = W-1.0; U<W+1.0; U+=0.05)
{
  cpa.SetParams( U, W, T );
  cpa.mu = U/2.0; //this is necessary for CalcDelta!!!!!!!!!!!!
  
  //set uniform mus between mu-W/2 and mu+W/2 
  cpa.DiscretizeMus(a,cpa.mu);
  
  //==========//

  //run IPT
  cpa.PrintIntermediate = 0;
  cpa.IntermediateCutoff = 10.0; //default 10.0
  sprintf(FN,"%s.W%.3f.T%.3f.U%.3f",prefix,W,T,U);
  cpa.intermediateName = FN;
  cpa.PrintDiffs = true;

  cpa.SetLoopOptions(1000, 1e-8);
  cpa.UseIPT = true;
 
  bool err = cpa.Run(&a);
  
  if (err) 
  { //reset Delta 
    for(int i=0; i<Nepsilons; i++)
    for(int j=0; j<N; j++)
      a.r[i].Delta[j] = 0.0;
    //go to the next W
    break;
  }

  sprintf(cmd,"mkdir %s",FN);
  system(cmd);
  a.PrintAllUberMinimal(FN,10.0);
  char fn[300];
  sprintf(fn,"%s/Delta",FN);
  PrintFunc(fn, N, a.r[0].Delta, a.r[0].omega);

  //printout average - should be put in IAResult
  complex<double>* G = new complex<double>[N];
  complex<double>* Sigma = new complex<double>[N];

  a.GetAverageGandSigma(G,Sigma);
  
  sprintf(fn,"%s/G_avg",FN);
  PrintFunc(fn, N, G, a.r[0].omega);
  sprintf(fn,"%s/Sigma_avg",FN);
  PrintFunc(fn, N, Sigma, a.r[0].omega);

  for(int i=0; i<N; i++)
        Sigma[i] = ii*a.r[0].omega[i] + U/2.0 - a.r[0].Delta[i] - 1.0/G[i];

  sprintf(fn,"%s/Sigma_eff",FN);
  PrintFunc(fn, N, Sigma, a.r[0].omega);

  delete [] G;
  delete [] Sigma;

  
/*
   //=================================//
  //-------prepare ctqmcSIAM-----//
  ctqmcSIAM siam;
  siam.M=20000000;  

  siam.runParallel = true;

  char rfn[300];
  sprintf(rfn,"/nfs/jaksa/QMC_statDMFT/met_qmc_runs/");
  siam.runFolderName = rfn;
  siam.execName = "/nfs/jaksa/QMC_statDMFT/ctqmc";

  cpa.ctqmcsiam = &siam;

  //run CTQMC starting from IPT
  omp_set_num_threads(1);

  cpa.SetPrintOutOptions(true, false);
  cpa.SetLoopOptions(20, 1e-5);
  cpa.UseIPT = false;

  sprintf(FN,"ctqmcCPA.W%.3f.T%.3f.U%.3f",W,T,U);
  cpa.intermediateName = FN;

  cpa.Run(&a);

  sprintf(cmd,"mkdir %s",FN);
  system(cmd);
  a.PrintAllMinimal(FN);
*/
  //print out mus 
 /* double* mus = new double[Nsites];
  dhm.GetMus(a, mus);
  char musFN[300];
  sprintf(musFN,"%smus",FN);
  PrintFunc(musFN,Nsites,mus);
  delete [] mus;
*/
  //global filling
/*  char gnFN[300];
  sprintf(gnFN,"%s/global_n",FN);
  FILE* gnFile = fopen(gnFN,"w");
  fprintf(gnFile,"%.15le\n",a.Global_n());
  fclose(gnFile);
*/
  //Pade on Green's function and Self-energy
/*  for(int id=0; id<Nsites; id++)
  { char pFN[300];
    sprintf(pFN,"%sGw.%d",FN,id);
    PadeToFile( 2000, a.r[id].G,  a.r[id].omega, pFN, 500, 3.0 );
    sprintf(pFN,"%sSigw.%d",FN,id);
    PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 500, 3.0 );

    printf("pade, site id=%d : DONE!\n",id);
  }
*/
}
}
}
  return 0;
}
