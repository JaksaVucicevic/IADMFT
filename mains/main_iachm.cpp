#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"

int main()
{

  // number of matsubara frequences to be used.
  // the fast fourier transform requires that this is a power of 2,
  // and that the number of tau points is the same
  int N = pow(2,13);

  // initialize the object that does the calculation
  IACHM iachm;

  // set DMFT loop option
  iachm.SetBroydenOptions(false, false, 0); 
  iachm.SetLoopOptions(100, 1e-9);

  // prepare the file where lambda will be output as a function of U and T  
  char lnUTFN[300] = "lambda.UT";
  FILE* lnUTFile = fopen(lnUTFN,"w");
  fclose(lnUTFile);

  // iterate over parameters
  for(double T=0.01; T<0.501; T+=0.01)
  {
    // the imag axis greid depends on temeprature. whenever temperature is changed, reinitialize grid
    IAGRID iagrid(N,N, T);
    // create the object for storing all the results. tau and omega arrays are automatically filled in by iagrid
    IAResult iaresult(&iagrid);
    // initialize FFT. this is also dependent on temperature
    iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau); 
   
    for(double n=0.5; n<0.7; n+=1000.05)
    {
      // if non interacting DOS is symetric around w=0 and n=0.5 then the calculation is simplified. use particle hole symmetry
      if (n!=0.5) iachm.PHSymmetricCase = false;
      else iachm.PHSymmetricCase = true;

      iaresult.n=n;

      for(double U=0.2; U<4.1; U+=0.025)
      { char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
        iachm.SetParams(U, T, 0.5);
        // lambdas and diffs vs iteration index will be output by the LambdaCalculator inside iachm
        iachm.LC->SetOutputFileName(bareFN);         
  
        //--------- RUN THE CODE---------//
        iachm.Run(&iaresult);
       
        // gather the last lambda and the best lambda.
        // the best is determined from the sequence of 5 iterations where it is most constant. this one should be used
        // the last is determined from the last 3 iterations
        double best_lambda = iachm.LC->best_lambda;
        double last_lambda = iachm.LC->lambdas[0];
        printf("best_lambda in main: %.3f\n",lambda);

        // printout lambda
        lnUTFile = fopen(lnUTFN,"a");    
        fprintf(lnUTFile,"%.15le %.15le %.15le %.15le %.15le\n", n, U, T, best_lambda, last_lambda);
        fclose(lnUTFile);

        // printout result
        char FN[300];
        sprintf(FN, "IACHM%s",bareFN);
        iaresult.PrintResult(FN);
        
        // use Pade to analytically continue Green's funtion to the real axis
        sprintf(FN, "Gw%s",bareFN);
        PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );
      }
      // add an empty line in the file so that pm3d can be used in gnuplot
      lnUTFile = fopen(lnUTFN,"a");     
      fprintf(lnUTFile,"\n");
      fclose(lnUTFile);
    }
  }
  fclose(lnUTFile);
}

