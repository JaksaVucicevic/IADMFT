#include "IASIAM.h"
#include "routines.h"
#include "Broyden.h"
#include "IAGRID.h"
#include "IAResult.h"
#include "Input.h"
#include "FFT.h"

#ifdef _OMP
#include <omp.h>
#endif

//================== Constructors/DEstructors ====================//

void IASIAM::Defaults()
{
  U = 2.0;
  T = 0.05;
  epsilon = 0;
  
  //broyden parameters
  MAX_ITS = 100; //default 100
  Accr = 1e-9; //default 1e-9

  // mu0 search
  max_tries = 2;
  UseBroydenFormu0 = true; 

  PatchTailWithAtomicLimit = true;
  AtomicCutoff = 12.0;

  AmoebaScanStart = -2.0;
  AmoebaScanEnd = 2.0; 
  AmoebaScanStep = 0.2;
  AmoebaMaxIts = 100;
  AmoebaForceScanAndPrintOut = false;
}

IASIAM::IASIAM()
{
  Defaults();
}

IASIAM::IASIAM(const char* ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- IASIAM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);

  input.ReadParam(U,"IASIAM::U");
  input.ReadParam(T,"IASIAM::T");
  input.ReadParam(epsilon,"IASIAM::epsilon");

  input.ReadParam(MAX_ITS,"IASIAM::MAX_ITS");
  input.ReadParam(Accr,"IASIAM::Accr");

  input.ReadParam(UseBroydenFormu0,"IASIAM::UseBroydenFormu0");
  input.ReadParam(max_tries,"IASIAM::max_tries");
  input.ReadParam(AmoebaScanStart,"IASIAM::AmoebaScanStart");
  input.ReadParam(AmoebaScanEnd,"IASIAM::AmoebaScanEnd");
  input.ReadParam(AmoebaScanStep,"IASIAM::AmoebaScanStep");
  input.ReadParam(AmoebaMaxIts,"IASIAM::AmoebaMaxIts");
  input.ReadParam(AmoebaForceScanAndPrintOut,"IASIAM::AmoebaForceScanAndPrintOut");
}

IASIAM::~IASIAM()
{

}

//========================= INITIALIZERS ===========================//

void IASIAM::SetT(double T)
{
  this->T = T;
}

void IASIAM::SetEpsilon(double epsilon)
{
  this->epsilon = epsilon;  
}

void IASIAM::SetU(double U)
{
  this->U = U;
}

void IASIAM::SetUTepsilon(double U, double T, double epsilon)
{
  this->T = T;
  this->epsilon = epsilon;  
  this->U = U;
}


void IASIAM::SetBroydenParameters(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

//========================= RUN IASIAM EITH FIXED Mu ==========================//

bool IASIAM::Run(IAResult* r) //output
{  
  this->r = r;
  N = r->iagrid->get_N(); 
  iagrid = r->iagrid;

  printf("    -------IASIAM: mu=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", r->mu, U, T, epsilon);
   
  complex<double>* V = new complex<double>[1];
  //----- initial guess ------// 
  V[0] = r->mu0; 
  //---------------------------//

  //------ SOLVE IASIAM ---------//
  double mu0inits [] = {0.0, 1.0, -1.0, -0.8, 2.0, 1.5, -1.5,  2.5, -2.5, 
                        -2.0, 0.05, 0.8, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5, -0.05, 
                        0.4, -0.4, 0.6, -0.6, 2.3, -2.3, 2.8, -2.8, 1.8, -1.8  }; 
  bool failed = !UseBroydenFormu0;  
  int c = 0;
  while ( UseBroydenFormu0 and ( UseBroyden<IASIAM>(1, MAX_ITS, Accr, &IASIAM::SolveSiam, this, V) != 1 ) )
  { c++;
    if ( (c >= max_tries) or (c >= sizeof(mu0inits)/sizeof(double) - 1 ) )
    {
      printf("\n\n\n\n==== IASIAM ERROR ====: Broyden mu0 search failed to converge. Now switching to amoeba ...\n\n\n\n");
      failed = true;  
      break;
    }
    V[0] = mu0inits[c]; 
    
    printf("==================== ====================== ========== TRYING new mu0 int!!! c = %d, mu0init = %f\n\n\n",c, mu0inits[c]);
  };
  if (failed) 
  {  V[0] = mu0inits[0];  
     Amoeba(Accr, V, &IASIAM::SolveSiam); 
  }
  delete [] V;
  //----------------------------//
  r->mu0 = mu0;

  printf("        IASIAM: mu0 = %f, n = %f\n",r->mu0, r->n);

  return false;
}

//========================== RUN IASIAM With FIXED n ==============================//
// applicable ONLY in solving Clean Hubbard Model which implies epsilon = 0 and NIDOS is needed on input.
//NOTE that MPT Bs will ALWAYS be one iteration late. They will converge to their real values
//when the DMFT loop converges. First DMFT Iteration is ALWAYS solved WITHOUT MPT Bs.

//TODO in case of asym NIDOS, mu and mu0 are not known EVEN FOR n=0.5 !!!!

bool IASIAM::Run_CHM(IAResult* r) //output
{  
  this->r = r;
  N = r->iagrid->get_N();
  iagrid = r->iagrid;

  epsilon = 0;
 
  printf("    ------- IASIAM for CHM: n=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", r->n, U, T, epsilon);
  
  //----------------- CALCULATION ----------------------//
  if (PHSymmetricCase)
  {
    mu0 = 0.0;
    get_G0();
  }
  else
  { 
    complex<double>* V = new complex<double>[1];
    //------initial guess---------//
    V[0] = r->mu0; //initial guess is always the last mu0. in first DMFT iteration it is 0
    //---------------------------//

    printf("         IASIAM: about to calc mu0. at the moment: mu0 = %.3f mu=%.3f\n",r->mu0, r->mu);
    double initGuesses [] = {-1.0, 1.0, 0.3, -0.3, 0.1, -0.1, 
                             -0.8, 0.8, -0.6, 0.6, -0.7, 0.7,
                             -3.0, 3.0, 0.9, -0.9, 0.05, -0.05, 
                              0.5, -0.5, 0.2, -0.2, 2.0, -2.0};
    
    bool converged = false; 

    int i;
    for (i=0; ( (i<sizeof(initGuesses)/sizeof(double)) and (i<max_tries) ); i++)
    {  printf("------ IASIAM: trying with init guess: %f\n",real(V[0]));
       converged = UseBroyden<IASIAM>(1, 50, 1e-8, &IASIAM::get_G0, this, V);  
     if (converged) break;
     else V[0] = initGuesses[i];
    }
    if ((i==max_tries)and(!converged))
    {  V[0] = initGuesses[0]; 
       Amoeba(Accr, V, &IASIAM::get_G0); 
    }
    delete [] V;
  }
  r->mu0 = mu0;
  //PrintFunc("G0", N, r->G0, r->omega);
  printf("    mu0 = %f\n", mu0);
  
  fft->FtoT(r->G0, r->G0_tau);
  get_SOCSigma_tau();
  fft->TtoF(r->SOCSigma_tau, r->SOCSigma);

  if (PHSymmetricCase)
  { r->mu = 0.5*U;
    r->n=0.5;
    get_G();
  }
  else
  { 
    complex<double>* V = new complex<double>[1];
    V[0] = r->mu;

    double initGuesses [] = {-1.0, 1.0, -0.8, 0.8, -0.6, 0.6, 0.2, -0.2, 2.0, -2.0};
    
    bool converged = false; 

    int i;
    for (i=0; ( (i<sizeof(initGuesses)/sizeof(double)) and (i<max_tries) ); i++)
    {  printf("------ IASIAM: trying with init guess: %f\n",real(V[0]));
       converged = UseBroyden<IASIAM>(1, MAX_ITS, 1e-8, &IASIAM::get_G, this, V);
       if (converged){ break; }
       else V[0] = initGuesses[i];
       
    }
    if ((i==max_tries)and(!converged))
    {  V[0] = initGuesses[0]; 
       Amoeba(Accr, V, &IASIAM::get_G); 
    }
    delete [] V;
  }
  printf("    mu = %f\n", r->mu);

  //-----------------------------------------------------//

  return false;
}

//=================================== FUNCTIONS ===================================//
/*
complex<double> IASIAM::AtomicLimitG(int i)
{
  return (1.0-r->n)/(ii*r->omega[i]+r->mu) + r->n/(ii*r->omega[i]+r->mu-U);
}

void IASIAM::PatchAtomicLimitSigma()
{
  for (int i=0; i<N; i++) 
    if (r->omega[i]>12.0) r->Sigma[i] = ii*r->omega[i] + r->mu - 1.0/AtomicLimitG(i);
}
*/

void IASIAM::PatchAtomicLimitSigma()
{
  for (int i=0; i<N; i++) 
    if (r->omega[i]>AtomicCutoff) r->Sigma[i] = AtomicLimitSigma(r->omega[i], r->mu, r->n, U);
}

void IASIAM::PatchAtomicLimitG()
{
  for (int i=0; i<N; i++) 
    if (r->omega[i]>AtomicCutoff) r->G[i] = AtomicLimitG(r->omega[i], r->mu, r->n, U);
}


void IASIAM::get_G0()
{
  for (int i=0; i<N; i++) 
    r->G0[i] = 1.0 / ( ii*r->omega[i] + mu0 - r->Delta[i] ); 
}

void IASIAM::get_G0(complex<double>* V)
{
  mu0 = real(V[0]);

  get_G0();
  
  r->n0 = get_n(r->G0);
  
  V[0] = mu0 + r->n0 - r->n;
  printf("get_G0: V[0] = %.5f\n",real(V[0]));
} 


double IASIAM::get_n(complex<double> X[])
{
  double sum = 0;
  for (int i=0; i<N; i++) 
    sum += real(X[i]);
  
  //!!!!!!!! 2 or not 2
  return 2.0*T*sum+0.5; 
}

void IASIAM::get_SOCSigma_tau()
{
  for (int m=0; m<N; m++)
   r->SOCSigma_tau[m] = sqr(U) * r->G0_tau[m] * r->G0_tau[m] * r->G0_tau[N-1-m];
}


double IASIAM::get_b()
{ //we used mu0 as (mu0 - epsilon - U*n) in G0, now we're correcting that  
  if (not PHSymmetricCase)
    return ( (1.0 - r->n)*U - r->mu + mu0 + epsilon ) 
           /    ( r->n * (1.0 - r->n) * sqr(U) );
  else 
    return 0.0;
}

void IASIAM::get_Sigma()
{
  if (not PHSymmetricCase)
  { printf("going through asymmetric\n");
    double b = get_b();    
    for (int i=0; i<N; i++) 
      r->Sigma[i] =  U*r->n + r->SOCSigma[i] 
                              / ( 1.0 - b * r->SOCSigma[i] );
    
  }
  else
  { printf("going through symmetric\n");
    for (int i=0; i<N; i++) 
      r->Sigma[i] =  U * r->n + r->SOCSigma[i];
  }

  if (PatchTailWithAtomicLimit) PatchAtomicLimitSigma();
}

//---------------- Get G -------------------------------//

void IASIAM::get_G()
{
  get_Sigma();

  for (int i=0; i<N; i++) 
    r->G[i] =  1.0 / (ii*r->omega[i] + r->mu - epsilon - r->Delta[i] - r->Sigma[i]) ;   

  if (PatchTailWithAtomicLimit) PatchAtomicLimitG();    
}

void IASIAM::get_G(complex<double>* V)
{
  r->mu = real(V[0]);

  get_G();

  V[0] = r->mu + get_n(r->G) - r->n;
} 

//---------------------------------------------------------------//

void IASIAM::SolveSiam(complex<double>* V)
{
  mu0 = real(V[0]);

  //--------------------//
  get_G0();

  r->n0 = get_n(r->G0); 
  r->n = r->n0;

  fft->FtoT(r->G0, r->G0_tau);
  get_SOCSigma_tau();
  fft->TtoF(r->SOCSigma_tau, r->SOCSigma);

  get_G();
   
  r->n = get_n(r->G);
  //--------------------//

  V[0] = mu0 + (r->n - r->n0); //we need to satisfy (get_n(G) == n) and 
}

void IASIAM::Amoeba(double accr, complex<double>* V, void (IASIAM::*f)(complex<double>*))
{
  //x here stands for mu0
  
  double x_start = AmoebaScanStart;
  double x_end   = AmoebaScanEnd;
  double x_step  = AmoebaScanStep;

  int sign_old=0;
  double x;
  bool found = false;
  int try_count = 0;
  double x_candidate;
  double x_best=0, diff_best=10e+100;
  
  while( (not found) and (try_count<1) ) 
  {
     FILE* ScanFile;  
     if (AmoebaForceScanAndPrintOut)
     {
       char ScanFN[50];
       sprintf(ScanFN, "scan.eps%.3f",epsilon);
       ScanFile = fopen(ScanFN,"w");
     }
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       (this->*f)(V);
       
       if (AmoebaForceScanAndPrintOut)
       {  char FN[50];
          sprintf(FN,"siam.eps%.3f.mu0_%.3f",epsilon, mu0);
          r->PrintResult(FN);
       }
      
       double x_res=real(V[0]);
       printf("         mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f r->n: %.3f\n",
          x, x-x_res, x_step, get_n(r->G0)- get_n(r->G),get_n(r->G),r->n);

       if (AmoebaForceScanAndPrintOut)
         fprintf(ScanFile,"%.15le %.15le %.15le %.15le\n", x, x-x_res, r->n, r->n0);

       if (sign_old==0) 
       { sign_old = int_sign(x-x_res);
         continue;
       }

       int sign = int_sign(x - x_res);
       if (abs(x-x_res) < diff_best) { x_best=x; diff_best = abs(x-x_res); };
       if ((sign_old!=sign) and (not found))
       {  x_candidate = x-x_step;
          found = true; 
          if (not AmoebaForceScanAndPrintOut) break; 
       }
    }
    try_count++;
    if (not found) { x_start *=2.0; x_end *= 2.0; x_step *= 2.0; printf("              mu0 candidate NOT found! now scanning a wider range...\n"); }
    if (AmoebaForceScanAndPrintOut) fclose(ScanFile);
  } 
 
  
  if (not found)
  {  printf("              mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
     V[0] = x_best;
     (this->*f)(V);
  }
  else
  {
    printf("              mu0 candidate found! proceeding with aomeba...\n");  
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    { it ++;
      V[0] = x;
      (this->*f)(V);
      double x_res=real(V[0]);
      printf("         it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(r->G0), get_n(r->G));
      converged = ( abs(x-x_res) < accr );
      int sign = int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) printf("          Amoeba: desired accuracy reached!\n");
  }
  printf("         --- Amoeba DONE ---\n");
}


