#include "../source/IAresArray.h"
#include "../source/IAResult.h"
#include "../source/IAGRID.h"
#include "../source/arrayInitializers.h"
#include "../source/IADHM.h"
#include "../source/routines.h"
#include "../source/pade.h"

#include <cstdlib>
#include <cstdio>

int main()
{
  //---- init non-interacting Hamiltonian ------//

  // square lattice 10x10
  int size = 5;
  int Nsites = sqr(size); 
  //int Nsites = pow(size,3); //cubic lattice
  double t=0.5; 

  double** H0 = Array2D< double >( Nsites, Nsites );
  initCubicTBH(size, size, 1, t, H0);
  //initCubicTBH(size, size, size, t, H0); //cubic lattice

  //PrintMatrix("H0", Nsites, Nsites, H0);
  
  //--------------- StatDMFT -------------------//
  //params
  double T=0.05;
  double U=6.0;
  double W=2.0;

  double mu=U/2.0; 
  
  //matsubara grid
  int Nw = pow(2,13);
  IAGRID g( Nw, Nw, T );

  //all the storage arrays
  IAresArray a( Nsites, &g );

  //-------------//
  IADHM iadhm;
  iadhm.H0 = H0;
  iadhm.SetParams( U, W, T );

  //init ffts
  iadhm.fft = new FFT[Nsites];
  for(int id=0; id<Nsites; id++)
    iadhm.fft[id].Initialize(Nw, T, a.r[id].omega, a.r[id].tau);   
  
  //set random mus between mu-W/2 and mu+W/2
  iadhm.RandomizeMus(a, mu, 1234); 
  //or set mus from an existing array
  //iadhm.SetMus(a, mus); 


  //printout mus in a file to check
  //for(int id=0; id<Nsites; id++)
  //  printf("id: %.3d mu: %.5f\n", id, a.r[id].mu);
  double* mus = new double[Nsites];
  iadhm.GetMus(a, mus);
  PrintFunc("mus",Nsites,mus);
  delete [] mus;
  
  //==========//
  iadhm.Run(&a);
  //==========//

  //print out results
  char FN[300];
  sprintf(FN,"IADHM.U%.3f.W%.3f.T%.3f/",U,W,T);
  char cmd[300];
  sprintf(cmd,"mkdir %s",FN);
  system(cmd);
  a.PrintAll(FN);

  //Pade on Green's function and Self-energy
/*  for(int id=0; id<Nsites; id++)
  { char pFN[300];
    sprintf(pFN,"%sGw.%d",FN,id);
    PadeToFile( 2000, a.r[id].G,  a.r[id].omega, pFN, 500, 3.0 );
    sprintf(pFN,"%sSigw.%d",FN,id);
    PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 500, 3.0 );
  }
*/
  //release H0
  FreeArray2D<double>(H0, Nsites);
  return 0;
}
