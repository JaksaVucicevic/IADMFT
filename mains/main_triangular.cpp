#include <cstdlib>
#include <cstdio>
#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"

int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double T = atof(argv[1]);

  int Nw = 1000;
  double* NIDOS = new double[Nw];
  double* w = new double[Nw];
  ReadFunc("triangular_dos", Nw, NIDOS, w);

  int N = pow(2,13);

  IACHM iachm;
  iachm.SetUseBethe(false);
  iachm.SetNIDOS(Nw, NIDOS, w);
  iachm.SetBroydenOptions(false, false, 0); 
  iachm.UseLambdaCalculator = true;
  iachm.SetLoopOptions(100, 1e-9);
  iachm.PHSymmetricCase = false;
/*  
  char lUTFN[300] = "lambda.UT";
  FILE* lUTFile = fopen(lUTFN,"w");
  fclose(lUTFile);

  FILE* instablFile = fopen("instabline","w");
  fclose(instablFile);
*/
  //for(double T=0.01; T<0.501; T+=0.01)
  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);
    iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau);

    double n=0.5;
    iaresult.n=n;

    char lTFN[300];
    sprintf(lTFN,"lambda.T%.3f",T);
    FILE* lTFile = fopen(lTFN,"w");
    fclose(lTFile);

    double min_best_lambda=1000;
    double max_best_lambda=-1000;
    double min_last_lambda=1000;
    double max_last_lambda=-1000;
    double Umin_best;
    double Umin_last;
    double Umax_best;
    double Umax_last;


    for(double U=0.2; U<5.0; U+=0.025)
    {   char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
        iachm.SetParams(U, T, 0.5);
        iachm.LC->SetOutputFileName(bareFN);         
  
        iachm.Run(&iaresult);

        double best_lambda = iachm.LC->best_lambda;
        double last_lambda = iachm.LC->lambdas[0];
        //printf("best_lambda in main: %.3f\n",lambda);
        
        /*lUTFile = fopen(lUTFN,"a");    
        fprintf(lUTFile,"%.15le %.15le %.15le %.15le\n", U, T, best_lambda, last_lambda);
        fclose(lUTFile);
        */
        lTFile = fopen(lTFN,"a");
        fprintf(lTFile,"%.15le %.15le %.15le\n", U, best_lambda, last_lambda);
        fclose(lTFile);
 
        if (best_lambda<min_best_lambda) { min_best_lambda = best_lambda; Umin_best = U; }      
        if (best_lambda>max_best_lambda) { max_best_lambda = best_lambda; Umax_best = U; }      
        if (last_lambda<min_last_lambda) { min_last_lambda = last_lambda; Umin_last = U; }      
        if (last_lambda>max_last_lambda) { max_last_lambda = last_lambda; Umax_last = U; }      
        
        char FN[300];
        sprintf(FN, "IACHM%s",bareFN);
        iaresult.PrintResult(FN);
        //sprintf(FN, "Gw%s",bareFN);
        //PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );
    }
    /*lUTFile = fopen(lUTFN,"a");     
    fprintf(lUTFile,"\n");
    fclose(lUTFile);    
    */

    char instablFN[300];
    sprintf(instablFN,"instabline.T%.3f",T);
    FILE* instablFile = fopen(instablFN,"w");
    fprintf(instablFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", T, Umin_best, Umax_best, Umin_last, Umax_last, min_last_lambda);
    fclose(instablFile);
    
  }
  //fclose(lUTFile);

  delete [] NIDOS;
  delete [] w;
}

