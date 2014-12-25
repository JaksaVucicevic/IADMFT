#include <cstdio>
#include "IACHM.h"
#include "IAGRID.h"
#include "IAResult.h"
#include "Input.h"
#include "IASIAM.h"
#include "routines.h"
#include <omp.h>

class IASIAM;

void IACHM::Defaults() 
{
    U = 2.0;
    T = 0.05;
    UseBethe = true;
    PHSymmetricCase = true;
    UseFixedMuSIAMRun = false;

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

  if (UseFixedMuSIAMRun)
    return iasiam->Run(r);
  else
    return iasiam->Run_CHM(r);
}

void IACHM::CalcDelta()
{
  if (UseBethe) printf("\n\n-- INFO -- CHM Self Consistency - Delta = t^G, t=%f\n\n",t);
  if (UseBethe)
    for (int i=0; i<N; i++)        
      r->Delta[i] = sqr(t) * r->G[i];        
  else    
    for (int i=0; i<N; i++)
    { complex<double>* g = new complex<double>[Nw];
      for(int j=0; j<Nw; j++)
        g[j] = NIDOS[j] / (ii*r->omega[i] + r->mu - w[j] - r->Sigma[i]);
      r->G[i] = TrapezIntegral(Nw, g, w);  
      r->Delta[i] = ii*r->omega[i] + r->mu - r->Sigma[i] - 1.0/r->G[i];    
      delete [] g;
    }

  if (PrintIntermediate)
  { char FN[300];
    sprintf( FN, "IACHM.U%.3f.T%.3f.it%d", U, T, Iteration);
    r->PrintResult(FN);
  }  
  
}
