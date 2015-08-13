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

int main(int argc, char* argv [])
{
  if(argc<4) exit(0);
  double W = atof(argv[1]);
  double T = atof(argv[2]);
  double U = atof(argv[3]);
  double mu=U/2.0; 

  //---- init non-interacting Hamiltonian ------//

  // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=0.5; 

  double** H0 = Array2D< double >( Nsites, Nsites );
  //initCubicTBH(size, size, 1, t, H0);
  initCubicTBH(size, size, size, t, H0); //cubic lattice

  //PrintMatrix("H0", Nsites, Nsites, H0);
  
  //--------------- StatDMFT -------------------//
  char cmd[300];   


  //matsubara grid
  IAGRID g2( 2000, 2000, T );

  //all the storage arrays
  IAresArray a2( Nsites, &g2 );

  //a2.ReadFromFiles("MiloshResult");

  for(int i=0; i<Nsites; i++)
  { 
    char fn[300];
    sprintf(fn,"MiloshResult/%d",i);
    int N,M;
    double** X;
    ReadFunc(fn, N, M, X);
    /*if (i%10==0)
    for(int k=1; k<M; k++)
    { char fn[300];
      sprintf(fn,"MiloshResultREAD/X_%d.%d",i,k);
      PrintFunc(fn, N, X[k], X[0]);
    }*/
    for(int j=0; j<N; j++)
    { //printf("X[0][%d]: %.3f, X[2][%d]: %.3f, X[3][%d]: %.3f, X[14][%d]: %.3f, X[15][%d]: %.3f\n", j, X[0][j],  j, X[2][j],   j, X[3][j],   j, X[15][j],   j, X[15][j] );
      a2.r[i].omega[j] = 3.0*X[0][j];
      a2.r[i].Delta[j] = complex<double>(X[2][j],X[3][j]);
      a2.r[i].G[j] = complex<double>(X[14][j],X[15][j]);
    }  
    for(int j=0; j<M; j++)
      delete [] X[j];
    delete [] X;
  }

  sprintf(cmd,"mkdir MiloshResultREAD");
  system(cmd);
  a2.PrintAll("MiloshResultREAD");
  //-------prepare ctqmcSIAM-----//
  ctqmcSIAM siam;
  siam.M=20000000;
  siam.runParallel = true;

  char rfn[300];
  sprintf(rfn,"/nfs/jaksa/QMC_statDMFT/");
  siam.runFolderName = rfn;
  siam.execName = "/nfs/jaksa/QMC_statDMFT/ctqmc";


  //-------------//
  ctqmcDHM dhm;
  dhm.ctqmcsiam = &siam;
  dhm.H0 = H0;
  dhm.SetParams( U, W, T );

  dhm.SetLoopOptions(20, 1e-4);
  dhm.SetBroydenOptions(false, false, 1e-4);

  //read epsilons from file
  double* epsilons = new double [Nsites];  
  ReadFunc("milosh.epsilons", Nsites, epsilons, epsilons);
  double* mus = new double [Nsites];    
  for(int i=0; i<Nsites; i++)
    mus[i] = mu - 3.0*epsilons[i];

  dhm.SetMus(a2,mus);
  
  //==========//
  omp_set_num_threads(1);

  dhm.SetPrintOutOptions(true, false);
  dhm.SetLoopOptions(20, 1e-5);
  dhm.UseIPT = false;

  //run CTQMC starting from milosh result
  char FN[300];
  sprintf(FN,"ctqmcDHM.fromMilosh.W%.3f.T%.3f.U%.3f",W,T,U);
  dhm.intermediateName = FN;

  dhm.Run(&a2);

  sprintf(cmd,"mkdir %s",FN);
  system(cmd);
  a2.PrintAllMinimal(FN);

  //==========//

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
  sprintf(gnFN,"%sglobal_n",FN);
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
  //release H0
  FreeArray2D<double>(H0, Nsites);

  return 0;
}
