#include "../source/IAresArray.h"
#include "../source/IAResult.h"
#include "../source/IAGRID.h"
#include "../source/arrayInitializers.h"
#include "../source/IADHM.h"
#include "../source/routines.h"
#include "../source/pade.h"

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
  int size = 2;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=0.5; 

  double** H0 = Array2D< double >( Nsites, Nsites );
  //initCubicTBH(size, size, 1, t, H0);
  initCubicTBH(size, size, size, t, H0); //cubic lattice

  PrintMatrix("H0", Nsites, Nsites, H0);
  
  //--------------- StatDMFT -------------------//
 
  //matsubara grid
  int Nw = pow(2,13);

  IAGRID g( Nw, Nw, T );

  //all the storage arrays
  IAresArray a( Nsites, &g );

  //-------------//
  IADHM iadhm;
  iadhm.H0 = H0;
  iadhm.SetParams( U, W, T );
  
  iadhm.PatchDelta = false;
  iadhm.SetLoopOptions(300, 1e-9);
  iadhm.SetBroydenOptions(false, false, 1e-4);

  //init ffts
  iadhm.fft = new FFT[Nsites];
  for(int id=0; id<Nsites; id++)
    iadhm.fft[id].Initialize(Nw, T, a.r[id].omega, a.r[id].tau);   

  //set random mus between mu-W/2 and mu+W/2
  iadhm.RandomizeMus(a, mu, 1234); 
   
  //==========//
  iadhm.Run(&a);
  //==========//

  //print out results
  char FN[300];
  sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f/",W,T,U);
  char cmd[300];
  sprintf(cmd,"mkdir %s",FN);
  system(cmd);
  a.PrintAll(FN);

  //print out mus 
  double* mus = new double[Nsites];
  iadhm.GetMus(a, mus);
  char musFN[300];
  sprintf(musFN,"%smus",FN);
  PrintFunc(musFN,Nsites,mus);
  delete [] mus;

  //Scattering rate histogram
  int Nbins = 12;
  double* X = new double[Nsites];

  for(int id=0; id<Nsites; id++)
  {  //X[id] = imag( pade( 2000, a.r[id].omega, a.r[id].Sigma, 0.0 ) );
     X[id] = imag( CubicFrom4points(a.r[id].Sigma,  a.r[id].omega) );
     printf("scattering rate cubic spline, site id=%d : DONE!\n",id);
  }

  double* x = new double[Nbins];
  double* P = new double[Nbins];
 
  Histogram(Nsites, X, Nbins, x, P);
  
  char thFN[300];
  sprintf(thFN,"%stau_histogram",FN);
  PrintFunc(thFN, Nbins, P, x);

  char tausFN[300];
  sprintf(tausFN,"%staus",FN);
  PrintFunc(tausFN, Nsites, X);

  delete [] x;
  delete [] P;
  delete [] X;

  //global filling
  char gnFN[300];
  sprintf(gnFN,"%sglobal_n",FN);
  FILE* gnFile = fopen(gnFN,"w");
  fprintf(gnFile,"%.15le\n",a.Global_n());
  fclose(gnFile);

  //Pade on Green's function and Self-energy
  /*#pragma omp parallel for
  for(int id=0; id<Nsites; id++)
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
