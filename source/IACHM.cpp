#include <cstdio>
#include "IACHM.h"
#include "IAGRID.h"
#include "IAResult.h"
#include "Input.h"
#include "IASIAM.h"
#include "routines.h"
#include "Mixer.h"

#include <omp.h>

class IASIAM;

void IACHM::Defaults() 
{
    U = 2.0;
    T = 0.05;
    UseBethe = true;
    PHSymmetricCase = true;
    UseFixedMuSIAMRun = false;

    PatchDelta = true;
    PatchTailWithAtomicLimit = false;

    t = 0.5;

    iasiam = new IASIAM();
}

IACHM::IACHM() : Loop()
{
  Defaults();
}

IACHM::IACHM(const char* ParamsFN) : Loop(ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- CHM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(U,"CHM::U");
  input.ReadParam(T,"CHM::T");
  input.ReadParam(t,"CHM::t");
  input.ReadParam(UseBethe,"CHM::UseBethe");
  input.ReadParam(UseFixedMuSIAMRun,"CHM::UseFixedMuSIAMRun");

  iasiam = new IASIAM(ParamsFN);
}

void IACHM::ReleaseMemory()
{
  printf("CHM release\n");
  iasiam->~IASIAM();
}

IACHM::~IACHM()
{
  ReleaseMemory();
}

void IACHM::SetParams(double U, double T, double t)
{
  this->U = U;
  this->T = T;
  this->t = t;
}

void IACHM::SetNIDOS(int Nw, double* NIDOS, double* w)
{
  this->Nw = Nw;
  this->NIDOS = NIDOS;
  this->w = w;
}

void IACHM::SetUseBethe(bool UseBethe)
{
  this->UseBethe = UseBethe;
}

//-------------------------------------------------------//
  
bool IACHM::SolveSIAM()
{
  //IASIAM iasiam;
  iasiam->SetUTepsilon(U,T,0.0);
  iasiam->PHSymmetricCase = PHSymmetricCase;
  iasiam->fft = &fft;
  
  iasiam->PatchTailWithAtomicLimit = PatchTailWithAtomicLimit;
  iasiam->AtomicCutoff = AtomicCutoff;

  if (UseFixedMuSIAMRun)
    return iasiam->Run(r);
  else
    return iasiam->Run_CHM(r);
}

//-----------------test iterative perturbation-------------------//
// experimental method.
// calculates the second order diagram iteratively o include a subclass of diagrams rather than a single diagram.

/*
#include "../source/pade.h"
bool IACHM::SolveSIAM()
{
  //IASIAM iasiam;
  iasiam->SetUTepsilon(U,T,0.0);
  iasiam->PHSymmetricCase = PHSymmetricCase;
  iasiam->fft = &fft;
  
  iasiam->PatchTailWithAtomicLimit = PatchTailWithAtomicLimit;
  iasiam->AtomicCutoff = AtomicCutoff;

  for(int n=0; n<iagrid->get_N(); n++)
    r->SOCSigma[n] = 0.0;    

  complex<double>* original_Delta = new complex<double>[r->iagrid->get_N()];
  for(int n=0; n<iagrid->get_N(); n++)
   original_Delta[n] = r->Delta[n];     

  bool err;
  int Nipt = 100;
  Mixer< complex<double> > mixer(iagrid->get_N(), 2, (int []) {1,0}, 1e-8);
  mixer.Mix(r->SOCSigma);

  for(int i=0; i<Nipt; i++)
  { printf("IPT iteration: %d\n",i);
     
    if (i!=0)
      for(int n=0; n<iagrid->get_N(); n++)
        r->Delta[n] = original_Delta[n] + r->SOCSigma[n];    

    if (UseFixedMuSIAMRun)
      err = iasiam->Run(r);
    else
      err = iasiam->Run_CHM(r);

    if (mixer.Mix(r->SOCSigma)) break; 

    char FN[300];
    sprintf(FN,"IPT.it%d",i);
    if (i%10==0) r->PrintResult(FN);
    //if(i==0)
    //  PadeToFile(2000, r->G,  r->omega, "Gw", 600, 4.0 );  
    
  }
  if (Nipt>1)
  for(int n=0; n<iagrid->get_N(); n++)
  { r->G[n] = 1.0/(ii*r->omega[n] + r->mu - original_Delta[n] - r->Sigma[n]);    
    //r->Sigma[n] = r->Delta[n] + r->Sigma[n];     
  }
  return err;
}
*/

void IACHM::CalcDelta()
{
  if (UseBethe) printf("\n\n-- INFO -- CHM Self Consistency - Delta = t^G, t=%f\n\n",t);
  if (UseBethe)
    for (int i=0; i<N; i++)        
      r->Delta[i] = sqr(t) * r->G[i];        
  else    
  { for (int i=0; i<N; i++)
    { complex<double>* g = new complex<double>[Nw];
      for(int j=0; j<Nw; j++)
        g[j] = NIDOS[j] / (ii*r->omega[i] + r->mu - w[j] - r->Sigma[i]);
      r->G[i] = TrapezIntegral(Nw, g, w);  
      delete [] g;
    }
    if (PatchTailWithAtomicLimit) r->PatchAtomicLimitG(AtomicCutoff, U);
    for (int i=0; i<N; i++)
      r->Delta[i] = ii*r->omega[i] + r->mu - r->Sigma[i] - 1.0/r->G[i];
    if (PatchDelta) r->PatchDelta(AtomicCutoff);
      
  }
  if (PrintIntermediate)
  { char FN[300];
    sprintf( FN, "IACHM.U%.3f.T%.3f.it%d", U, T, Iteration);
    r->PrintResult(FN);
  }  
  
}
