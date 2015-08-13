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


//============================================== Pade G and Delta, get Sigma from RASC, calc local resistivity and histogram================================================//


int main(int argc, char* argv [])
{
//  if(argc<4) exit(0);
//  double W = atof(argv[1]);
//  double T = atof(argv[2]);
//  double U = atof(argv[3]);
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

  double T=0.01;
  //for(double T=0.002; T<0.52; T*=2.0)
  { IAGRID g( Nw, Nw, T );

    //all the storage arrays
    IAresArray a( Nsites, &g );
    
    double W=6.0;
    //for(double W=5.0; W<7.1; W+=100.0)
    { 
      
    for(double U=1.5; U<7.6; U+=2.0)
    { //==========//

      int Nseeds = 195;

      double** mus = new double*[Nseeds];
      double** ns = new double*[Nseeds];
      double** Aw0s = new double*[Nseeds];
      double** Aw0s_pade = new double*[Nseeds];
      double** Zs_RASC = new double*[Nseeds];
      double** taus_RASC = new double*[Nseeds];

      for(int seed=0; seed<Nseeds; seed++)
      {
        //load the result and clear previous sorted files
        sprintf(FN,"IADHM.seed%d.W%.3f.T%.3f.U%.3f",seed,W,T,U);
      
        bool successful = a.ReadFromMinimalFiles(FN);
        if (not successful) return -1;
        printf("======================= FILE LOADED: seed:%d U: %.3f T: %.3f W: %.3f =======================\n",seed, U,T,W);           
      
        mus[seed] = a.Get_mus();
        //normalize mus to go from -1 to 1
        // for(int i=0; i<Nsites; i++)
        //{  mus[seed][i] -= U/2.0;
        //   mus[seed][i] /= W/2.0;
        //}
        ns[seed] = a.Get_ns();
        Aw0s[seed] = a.Get_Aw0s();
        Aw0s_pade[seed] = a.Get_pade_Aw0s();
        Zs_RASC[seed] = a.Get_Zs_RASC();
        taus_RASC[seed] = a.Get_taus_RASC();
      }
     
      int NXbins = 100;
      int NYbins = 15;
      double* x = new double[NXbins];
      double** y = Array2D<double>(NXbins, NYbins);
      double** P = Array2D<double>(NXbins, NYbins);    

      Histogram2D(Nseeds, Nsites, mus, ns, NXbins, NYbins, x, y, P, false);
      sprintf(pFN,"n_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, Aw0s, NXbins, NYbins, x, y, P, false);
      sprintf(pFN,"Aw0_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, Aw0s, NXbins, NYbins, x, y, P, true);
      sprintf(pFN,"Aw0_log_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, Aw0s_pade, NXbins, NYbins, x, y, P, false);
      sprintf(pFN,"Aw0_pade_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, Aw0s_pade, NXbins, NYbins, x, y, P, true);
      sprintf(pFN,"Aw0_pade_log_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, Zs_RASC, NXbins, NYbins, x, y, P, false);
      sprintf(pFN,"Z_RASC_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      Histogram2D(Nseeds, Nsites, mus, taus_RASC, NXbins, NYbins, x, y, P, true);
      sprintf(pFN,"tau_RASC_log_2Dhistogram.W%.3f.T%.3f.U%.3f",W,T,U);
      Print2DHistogram(pFN, NXbins, NYbins, x, y, P);

      FreeArray2D<double>(P, Nseeds);
      FreeArray2D<double>(y, Nseeds);
      delete [] x;

      FreeArray2D<double>(mus, Nseeds);
      FreeArray2D<double>(ns, Nseeds);
      FreeArray2D<double>(Aw0s, Nseeds);
      FreeArray2D<double>(Aw0s_pade, Nseeds);
      FreeArray2D<double>(Zs_RASC, Nseeds);
      FreeArray2D<double>(taus_RASC, Nseeds);
    }
  }
  }

  return 0;
}


