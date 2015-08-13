#include <cstdlib>
#include <cstdio>
#include "../source/ctqmcCHM.h"
#include "../source/ctqmcSIAM.h"
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

  //prepare ctqmcSIAM
  ctqmcSIAM siam;
  siam.M=100000000;
  siam.runParallel = false;

  char rfn[300];
  //sprintf(rfn,"/nfs/jaksa/QMC_triang/RUN%d",Nproc);
  sprintf(rfn,"/nfs/jaksa/QMC_triang/run_ctqmc.T%.3f",T);
  siam.runFolderName = rfn;
  siam.execName = "/nfs/jaksa/QMC_triang/ctqmc";

  ctqmcCHM chm;
  chm.ctqmcsiam = &siam;  
  chm.SetUseBethe(false);
  chm.SetNIDOS(Nw, NIDOS, w);
  chm.SetBroydenOptions(false, false, 0); 
  chm.UseLambdaCalculator = true;
  chm.PHSymmetricCase = false;
  chm.PrintMuHistory = true;
  chm.PatchDelta = true;
  chm.PatchTailWithAtomicLimit = false;
  chm.AtomicCutoff = 8.0;
  chm.SetPrintOutOptions(false, false);
  chm.UseBroydenForMu = false;


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
    chm.fft.Initialize(N, T, iaresult.omega, iaresult.tau);

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

    char muFN[300];
    sprintf(muFN,"mu_vs_U.T%.3f",T);
    FILE* muFile = fopen(muFN,"w");
    fclose(muFile); 

    double n = 0.5;
    for(double U=1.5; U<3.5; U+=0.1)
    {   char bareFN[300];
        sprintf(bareFN, ".n%.3f.U%.3f.T%.3f",n,U,T);
        chm.SetParams(U, T, 0.5);
        chm.LC->SetOutputFileName(bareFN);         
  
        //------- clear Delta ---------//
        for(int i=0; i<N; i++)
          iaresult.Delta[i] = 0;

        //------- IPT ---------//
        chm.SetLoopOptions(100, 1e-9);     

        chm.UseFixedMuSIAMRun = false;
        iaresult.n = n;
        iaresult.mu = U/2.1;   

        chm.LC->SetDoOutput(false);
        printf("about to run IPT\n");

        chm.UseIPT = true;
        chm.Run(&iaresult);
        char FN[300];
        sprintf(FN, "IACHM%s",bareFN);
        iaresult.PrintResult(FN);   

        double IPT_mu = iaresult.mu;

        //------ find mu ---------//
        chm.SetLoopOptions(15, 1e-5);     

        chm.UseIPT = false;
        printf("about to run CTQMC\n");
        chm.Run(&iaresult);
       
        muFile = fopen(muFN,"a");
        fprintf(muFile,"%.15le %.15le %.15le %.15le\n", U, iaresult.mu, iaresult.n0, IPT_mu);
        fclose(muFile);

        //------- clear Delta ---------//
        for(int i=0; i<N; i++)
          iaresult.Delta[i] = 0;
        
        //------ run with fixed mu to find lambda -------//
        chm.UseFixedMuSIAMRun = true;

        chm.LC->SetDoOutput(true);
        chm.Run(&iaresult);

        //----- printout stuff -------//
        double best_lambda = chm.LC->best_lambda;
        double last_lambda = chm.LC->lambdas[0];
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
        
        sprintf(FN, "ctqmcCHM%s",bareFN);
        iaresult.PrintMinimal(FN);
        //sprintf(FN, "Gw%s",bareFN);
        //PadeToFile(2000, iaresult.G,  iaresult.omega, FN, 600, 4.0 );

        char cmd[300];
        sprintf(cmd, "mv mu_history mu_history%s",bareFN);
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

