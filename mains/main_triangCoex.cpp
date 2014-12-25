#include "../source/IACHM.h"
#include "../source/IAGRID.h"
#include "../source/IAResult.h"
#include "../source/pade.h"
#include "../source/LambdaCalculator.h"
#include "../source/routines.h"
#include <limits>

int main()
{
  int Nw = 1000;
  double* NIDOS = new double[Nw];
  double* w = new double[Nw];
  ReadFunc("triangular_dos", Nw, NIDOS, w);

  int N = pow(2,13);

  IACHM iachm;
  iachm.SetUseBethe(false);
  iachm.SetNIDOS(Nw, NIDOS, w);
  iachm.SetBroydenOptions(false, false, 0); 
  iachm.SetLoopOptions(100, 1e-7);
  iachm.PHSymmetricCase = false;

  char sFN[300] = "spinodals";
  FILE* sFile = fopen(sFN,"w");
  fclose(sFile);
  for(double T=0.006; T>0.001; T-=0.002)
  { IAGRID iagrid(N,N, T);
    IAResult iaresult(&iagrid);
    iachm.fft.Initialize(N, T, iaresult.omega, iaresult.tau);
    for(int n=0; n<N; n++) iaresult.Delta[n] = 0.25*(0.5/(ii*iaresult.omega[n]+5.0) + 0.5/(ii*iaresult.omega[n]-5.0));
    iaresult.n=0.5;

    double Um = std::numeric_limits<double>::quiet_NaN();
    double Ui = std::numeric_limits<double>::quiet_NaN();

    double U_start = 10.0;
    double U_step_max = 0.5;
    double U_step_min = 0.001;
    double U_max = 15.0;
    double U_min = 1.0;
  
    double U_step = -U_step_max;
    bool found = false;
    for(double U=U_start; (U_step<0) ? (U>U_min) : (U<U_max); U+=U_step)
    { 
      IAResult iarcpy(iaresult);

      iachm.SetParams(U, T, 0.5);
      iachm.LC->SetDoOutput(false);         
  
      iachm.Run(&iaresult);

      char FN[300];
      sprintf(FN, "%s.T%.3f.U%.3f",(U_step>0)?"Met":"Ins",T,U);
      iaresult.PrintResult(FN);

       
      if ( (U_step>0)
           and
           ( abs( imag(iaresult.G[0])/imag(iarcpy.G[0]) ) < 0.5
           )
         )
      { 
        found = true;
        U -= U_step;
        iaresult.CopyFrom(iarcpy);
      }

      if ( (U_step<0)
           and
           ( abs( imag(iaresult.Sigma[0])/imag(iarcpy.Sigma[0]) ) < 0.5
           )
         )
      { 
        found = true;
        U -= U_step;
        iaresult.CopyFrom(iarcpy);

      }

      if (found) U_step /= 2.0;
      
      if (abs(U_step)<U_step_min)
      {
        if (U_step>0)
        { Um = U;
          break;
        }
        else
        { Ui = U;
          found = false;
          U-=2.0*U_step_max;
          U_step = U_step_max;
        } 
      }

    }
    sFile = fopen(sFN,"a");
    fprintf(sFile,"%.15le %.15le %.15le\n",T, Um, Ui);
    fclose(sFile);    
  }

} 

