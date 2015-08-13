#include "../source/arrayInitializers.h"
#include "../source/routines.h"
#include "../source/resistorNetwork.h"
#include <cstdlib>
#include <cstdio>

int main()
{
  // square lattice 10x10
  long long int size = 6;
  //int Nsites = sqr(size); 
  int Nsites = pow(size,3); //cubic lattice
 
  char FN[300];

  for(double W=2.5; W<7.6; W+=2.5)   
  {
    char RWFN[300];
    sprintf(RWFN,"R_vs_UT.W%.3f",W);
    FILE* RWFile = fopen(RWFN,"w");

    for(double U=W-1.0; U<W+1.0; U+=0.05)
    {
      char RWUFN[300];
      sprintf(RWUFN,"R_vs_T.W%.3f.U%.3f",W,U);
      FILE* RWUFile = fopen(RWUFN,"w");

      for(double T=0.002; T<0.52; T*=2.0)
      {
        sprintf(FN,"IADHM.W%.3f.T%.3f.U%.3f/rhos",W,T,U);
        printf(">>>>> Working %s\n",FN);
        double* rhos = new double[Nsites];
        double* mus = new double[Nsites];
     
        ReadFunc(FN, Nsites, rhos, mus);

        //for(int i=0; i<Nsites; i++)
        //  printf("i: %d rho: %.3le\n",i,rhos[i]);

        double Ravg = Resistivity(3, size, size, size, 1, 1, 1, rhos, false);
        double Rx = Resistivity(0, size, size, size, 1, 1, 1, rhos, false);
        double Ry = Resistivity(1, size, size, size, 1, 1, 1, rhos, false);
        double Rz = Resistivity(2, size, size, size, 1, 1, 1, rhos, false);
        double rho_typ = typical(Nsites, rhos);
        double rho_avg = average(Nsites, rhos);

        fprintf(RWFile,"%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n",U,T, Ravg, Rx, Ry, Rz, rho_typ, rho_avg);
        fprintf(RWUFile,"%.15le %.15le %.15le %.15le %.15le %.15le %.15le\n",T, Ravg, Rx, Ry, Rz, rho_typ, rho_avg);

        delete [] rhos;
        delete [] mus;
      }
      fprintf(RWFile,"\n");
      fclose(RWUFile);
    }
    fclose(RWFile);
  }
}
