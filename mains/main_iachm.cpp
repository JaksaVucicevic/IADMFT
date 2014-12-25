#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"

int main()
{
  int N = pow(2,13);

  IACHM iachm;

  iachm.SetBroydenOptions(false, false, 0); 
  iachm.SetLoopOptions(100, 1e-9);
  
  char lnUTFN[300] = "lambda.UT";
  FILE* lnUTFile = fopen(lnUTFN,"w");
  fclose(lnUTFile);
  for(double T=0.01; T<0.501; T+=0.01)
  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);
    iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau);
    for(double n=0.5; n<0.7; n+=1000.05)
    {
      if (n!=0.5) iachm.PHSymmetricCase = false;
      else iachm.PHSymmetricCase = true;
      iaresult.n=n;

      for(double U=0.2; U<4.1; U+=0.025)
      { char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
        iachm.SetParams(U, T, 0.5);
        iachm.LC->SetOutputFileName(bareFN);         
  
        iachm.Run(&iaresult);
        double best_lambda = iachm.LC->best_lambda;
        double last_lambda = iachm.LC->lambdas[0];
        printf("best_lambda in main: %.3f\n",lambda);

        lnUTFile = fopen(lnUTFN,"a");    
        fprintf(lnUTFile,"%.15le %.15le %.15le %.15le %.15le\n", n, U, T, best_lambda, last_lambda);
        fclose(lnUTFile);
        //char FN[300];
        //sprintf(FN, "IACHM%s",bareFN);
        //iaresult.PrintResult(FN);
        //sprintf(FN, "Gw%s",bareFN);
        //PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );
      }
      lnUTFile = fopen(lnUTFN,"a");     
      fprintf(lnUTFile,"\n");
      fclose(lnUTFile);
    }
  }
  fclose(lnUTFile);
}

