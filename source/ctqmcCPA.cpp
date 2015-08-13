#include "ctqmcCPA.h"
#include "ctqmcSIAM.h"
#include "IAResult.h"
#include "IAresArray.h"
#include "routines.h"

#include <cstdlib>
#include <cstdio>

ctqmcCPA::ctqmcCPA(): ctqmcDHM()
{
  UseBethe = true;
  UseAverageSigma = false;
  PHSymmetry = true;
  t = 0.5;
}


ctqmcCPA::~ctqmcCPA()
{ //was eva
} 

void ctqmcCPA::DiscretizeMus(IAresArray &a, double mu)
{
  double dW = ( (PHSymmetry) ? 0.5 : 1.0 ) 
              * W
              / ( (double) ( a.get_N()
                             - ( (PHSymmetry)?0:1 ) 
                           )
                );
  for(int i=0; i<a.get_N(); i++)
    a.r[i].mu = mu - W/2.0 + i*dW;  
}


void ctqmcCPA::CalcDelta()
{
  //mu here is set from outside. see ctqmcCPA.h   

  complex<double>* G = new complex<double>[N];
  complex<double>* Sigma = new complex<double>[N];  
  
  printf("-------CPA: about to average G and Sigma...\n"); 

  a->GetAverageGandSigma(G,Sigma);  

  printf("...done! Now about to calculate Delta!\n"); 

  #pragma omp parallel for
  for(int i=0; i<N; i++) 
  { if (PHSymmetry)
    { G[i]=complex<double>(0.0,imag(G[i]));
      Sigma[i]=complex<double>(U/2.0,imag(Sigma[i]));    
    }
 
    if (UseBethe)
      for(int j=0; j<Nsites; j++)
        a->r[j].Delta[i] = sqr(t)*G[i];
    else
    { 
      if (not UseAverageSigma)
        Sigma[i] = ii*a->r[0].omega[i] + mu - a->r[0].Delta[i] - 1.0/G[i]; 

      complex<double>* g = new complex<double>[Nw];
      for(int l=0; l<Nw; l++)
          g[l] = NIDOS[l] / (ii*a->r[0].omega[i] + mu - w[l] - Sigma[i]); 
      G[i] = TrapezIntegral(Nw, g, w);  //G becomes G^latt 
      delete [] g;
    
      if ((PatchDelta)and(i>0))
        if (a->r[0].omega[i-1]>AtomicCutoff) continue; // i-1 because we need one extra point. the delta patch starts from the first omega[i]>atomic cutoff.
      
      a->r[0].Delta[i] = ii*a->r[0].omega[i] + mu - Sigma[i] - 1.0/G[i];
      for(int j=1; j<Nsites; j++)
        a->r[j].Delta[i] = a->r[0].Delta[i];
    }
       
  } 
  printf("...done!\n"); 
  if (PatchDelta)
  { printf("---------CPA: Now about to patch Delta!\n"); 
    #pragma omp parallel for
    for(int j=0; j<Nsites; j++)
       a->r[j].PatchDelta(AtomicCutoff);
  }

  delete [] G;
  delete [] Sigma;
}
