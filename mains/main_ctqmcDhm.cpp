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
  if(argc<2) exit(0);
//  double W = atof(argv[1]);
//  double T = atof(argv[1]);
//  double U = atof(argv[3]);
  int seed = atoi(argv[1]);
  

  //---- init non-interacting Hamiltonian ------//

  // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=1.0/6.0; 

  double** H0 = Array2D< double >( Nsites, Nsites );
  //initCubicTBH(size, size, 1, t, H0);
  initCubicTBH(size, size, size, t, H0); //cubic lattice

  //PrintMatrix("H0", Nsites, Nsites, H0);
  
  //--------------- StatDMFT -------------------//
  ctqmcDHM dhm;
  dhm.H0 = H0;

  char FN[300];
  char cmd[300];


  //matsubara grid
  int Nw = pow(2,13);

  double T=0.01;
  //for(double T=0.002; T<0.52; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    //init ffts
    dhm.fft = new FFT[Nsites];
    for(int id=0; id<Nsites; id++)
      dhm.fft[id].Initialize(Nw, T, a.r[id].omega, a.r[id].tau);   

    //int seed = 1234;
    //for(int seed = 0; seed<100; seed+=1)
    {  //initialized Delta to (wide band) metallic
       
    double W=6.0;
    //for(double W=2.5; W<7.6; W+=2.5)
    { 
       for(int i=0; i<Nsites; i++)
       for(int j=0; j<Nw; j++)
         a.r[i].Delta[j] = sqr(7.0)*SemiCircleG(ii*a.r[i].omega[j],7.0);   
      
    for(double U=1.5; U<7.6; U+=2.0)
    //for(double U = W-1.0; U<W+1.0; U+=0.05)
    { double mu=U/2.0; 

      //set random mus between mu-W/2 and mu+W/2
      dhm.RandomizeMus(a, mu, seed); 

      dhm.SetParams( U, W, T );
  
      //==========//

      //run IPT
      dhm.PrintIntermediate = 30;
      dhm.IntermediateCutoff = 10.0; //default 10.0
      sprintf(FN,"IADHM.seed%d.W%.3f.T%.3f.U%.3f",seed,W,T,U);
      //sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      dhm.intermediateName = FN;
      dhm.PrintDiffs = true; 

      dhm.SetLoopOptions(1000, 1e-7);
      dhm.SetBroydenOptions(false,false,0);
      dhm.UseIPT = true;

      dhm.Run(&a);

      sprintf(cmd,"mkdir %s",FN);
      system(cmd);
      a.PrintAllMinimal(FN,10.0);     
    }
    }
    }
  }
  FreeArray2D<double>(H0, Nsites);
  return 0;
}


/*
int main(int argc, char* argv [])
{
//  if(argc<4) exit(0);
//  double W = atof(argv[1]);
//  double T = atof(argv[2]);
//  double U = atof(argv[3]);
//  double mu=U/2.0; 


  //---- init non-interacting Hamiltonian ------//

  // square lattice 10x10
  int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
  double t=1.0/6.0; 

  double** H0 = Array2D< double >( Nsites, Nsites );
  //initCubicTBH(size, size, 1, t, H0);
  initCubicTBH(size, size, size, t, H0); //cubic lattice

  //PrintMatrix("H0", Nsites, Nsites, H0);
  
  //--------------- StatDMFT -------------------//
  ctqmcDHM dhm;
  dhm.H0 = H0;

  char FN[300];
  char cmd[300];


  //matsubara grid
  int Nw = pow(2,13);

  for(double T=0.01; T<0.011; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );

    //init ffts
    dhm.fft = new FFT[Nsites];
    for(int id=0; id<Nsites; id++)
      dhm.fft[id].Initialize(Nw, T, a.r[id].omega, a.r[id].tau);   


    //initialized Delta to (wide band) metallic
    for(int i=0; i<Nsites; i++)
    for(int j=0; j<Nw; j++)
      a.r[i].Delta[j] = sqr(7.0)*SemiCircleG(ii*a.r[i].omega[j],7.0);

    //-------prepare ctqmcSIAM-----//
     //unnecessary if IPT
//    ctqmcSIAM siam;
//    siam.M=20000000;  

//    siam.runParallel = true;

//    char rfn[300];
//    sprintf(rfn,"/nfs/jaksa/QMC_statDMFT/met_qmc_runs/");
//    siam.runFolderName = rfn;
//    siam.execName = "/nfs/jaksa/QMC_statDMFT/ctqmc";

//    dhm.ctqmcsiam = &siam;   
    //-------------//

    
    for(double W=2.0; W<7.1; W+=1.0)
    { if ((W>4.8)and(W<5.2)) continue;
      char AwFN[300];
      sprintf(AwFN, "Aw0_vs_U.W%.3f.T%.3f",W,T);
      FILE* AwFile = fopen(AwFN,"w");
      fclose(AwFile); 
    for(double U=0.6; U<W+2.01; U+=0.1)
    { double mu=U/2.0; 
      dhm.SetParams( U, W, T );

      //set random mus between mu-W/2 and mu+W/2
      dhm.RandomizeMus(a, mu, 1234); 
  
      //==========//

      //run IPT
      dhm.PrintIntermediate = 30;
      dhm.IntermediateCutoff = 10.0; //default 10.0
      sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f",W,T,U);
      dhm.intermediateName = FN;
      dhm.PrintDiffs = true; 

      dhm.SetLoopOptions(1000, 1e-7);
      dhm.SetBroydenOptions(false,false,0);
      dhm.UseIPT = true;

      dhm.Run(&a);

      sprintf(cmd,"mkdir %s",FN);
      system(cmd);
      a.PrintAll(FN);
      a.PrintSortedMinimal(FN, 10.0);


      //A(w=0) histogram
      int Nbins = 12;
      double* X = new double[Nsites];

      for(int id=0; id<Nsites; id++)
      {  //X[id] = imag( pade( 2000, a.r[id].omega, a.r[id].Sigma, 0.0 ) );
        X[id] = -(1.0/pi)*imag( CubicFrom4points(a.r[id].G,  a.r[id].omega) );
        printf("A(w=0) cubic spline, site id=%d : DONE!\n",id);
      }

      double* x = new double[Nbins];
      double* P = new double[Nbins];
 
      Histogram(Nsites, X, Nbins, x, P);
  
      char thFN[300];
      sprintf(thFN,"%s/Aw0_histogram",FN);
      PrintFunc(thFN, Nbins, P, x);

      char tausFN[300];
      sprintf(tausFN,"%s/Aw0s",FN);
      PrintFunc(tausFN, Nsites, X);

      AwFile = fopen(AwFN,"a");
      for(int b=0; b<Nbins; b++)
        fprintf(AwFile,"%.15le %.15le %.15le\n",U,x[b],P[b]);
      fprintf(AwFile,"\n");
      fclose(AwFile); 


      //run CTQMC starting from IPT
//      /*  omp_set_num_threads(1); 

      dhm.SetPrintOutOptions(true, false);
      dhm.SetLoopOptions(20, 1e-5);
      dhm.UseIPT = false;

      sprintf(FN,"ctqmcDHM.fromMet.W%.3f.T%.3f.U%.3f",W,T,U);
      dhm.intermediateName = FN;

      dhm.Run(&a);

      sprintf(cmd,"mkdir %s",FN);
      system(cmd);
      a.PrintAllMinimal(FN);
//  
      //Pade on Green's function and Self-energy
//  /*for(int id=0; id<Nsites; id++)
      { char pFN[300];
        sprintf(pFN,"%sGw.%d",FN,id);
        PadeToFile( 2000, a.r[id].G,  a.r[id].omega, pFN, 500, 3.0 );
        sprintf(pFN,"%sSigw.%d",FN,id);
        PadeToFile( 2000, a.r[id].Sigma,  a.r[id].omega, pFN, 500, 3.0 );

        printf("pade, site id=%d : DONE!\n",id);
      }
//      
      //release H0
      
    }
  }
  }
  FreeArray2D<double>(H0, Nsites);
  return 0;
}
*/
