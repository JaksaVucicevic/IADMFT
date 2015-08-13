#include <cstdio>
#include <cstdlib>
#include "IADHM.h"
#include "IAGRID.h"
#include "IAResult.h"
#include "IAresArray.h"
#include "Input.h"
#include "IASIAM.h"
#include "routines.h"
#include "arrayInitializers.h"
#include "mkl_invertMatrix.h"

#ifdef _OMP
#include <omp.h>
#endif

class IASIAM;

void IADHM::Defaults() 
{
    U = 2.0;
    W = 1.0; 
    T = 0.05;

    //UseBethe = true;
    PHSymmetricCase = false;
    PatchDelta = false;
    PatchTailWithAtomicLimit = false;
    AtomicCutoff = 10.0;
}

IADHM::IADHM() : Loop()
{
  Defaults();
}

IADHM::IADHM(const char* ParamsFN) : Loop(ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- CHM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(U,"IADHM::U");
  input.ReadParam(T,"IADHM::T");
  input.ReadParam(W,"IADHM::W");
//  input.ReadParam(UseBethe,"IACHM::UseBethe");
//  input.ReadParam(UseFixedMuSIAMRun,"IACHM::UseFixedMuSIAMRun");

}

void IADHM::ReleaseMemory()
{
  printf("DHM release start\n");
  if (fft!=NULL) delete [] fft;
  printf("DHM release end\n");
}

IADHM::~IADHM()
{
  ReleaseMemory();
}

void IADHM::SetParams(double U, double W, double T)
{
  this->U = U;
  this->T = T;
  this->W = W;
}
/*
void IADHM::SetNIDOS(int Nw, double* NIDOS, double* w)
{
  this->Nw = Nw;
  this->NIDOS = NIDOS;
  this->w = w;
}

void IADHM::SetUseBethe(bool UseBethe)
{
  this->UseBethe = UseBethe;
}
*/
//-------------------------------------------------------//
  
void IADHM::RandomizeMus( IAresArray &a, double mu, unsigned int seed )
{
  srand (seed);
  int Nsites = a.get_N();
  
  double sum = 0;
  for(int id=0; id<Nsites; id++)
  {  a.r[id].mu = mu - W/2.0 + W * ((double) rand()) / ((double) RAND_MAX);
     sum += a.r[id].mu;
  }
  double avg_mu = sum/((double)Nsites);
  for(int id=0; id<Nsites; id++)
    a.r[id].mu += mu - avg_mu;
  
}

void IADHM::GetMus( IAresArray &a, double* mus )
{
  int Nsites = a.get_N();
  for(int id=0; id<Nsites; id++)
    mus[id] = a.r[id].mu;
}

void IADHM::SetMus( IAresArray &a, double* mus )
{
  int Nsites = a.get_N();
  for(int id=0; id<Nsites; id++)
    a.r[id].mu = mus[id];
}



bool IADHM::SolveSIAM()
{
  iasiams = new IASIAM[Nsites];
  #pragma omp parallel for
  for(int id=0; id<Nsites; id++)
  { iasiams[id].PatchTailWithAtomicLimit = PatchTailWithAtomicLimit;
    iasiams[id].AtomicCutoff = AtomicCutoff;
    iasiams[id].SetUTepsilon(U,T,0.0);
    iasiams[id].PHSymmetricCase = PHSymmetricCase;
    iasiams[id].fft = &(fft[id]);
    iasiams[id].Run(&(a->r[id]));
  }
  delete [] iasiams;
  return false;
}

void IADHM::CalcDelta()
{
  #pragma omp parallel for 
  for(int n=0; n<N; n++)
  { 
    //printf("IADMFT::CalcDelta: working iw_%d\n",n);
    complex<double>** invG = Array2D< complex<double> >(Nsites, Nsites);

    for(int i=0; i<Nsites; i++)
    for(int j=0; j<Nsites; j++)
      invG[i][j] = H0[i][j];

    for(int id=0; id<Nsites; id++)
      invG[id][id] += ii*a->r[id].omega[n] + a->r[id].mu - a->r[id].Sigma[n];   

    //if(n==0) PrintMatrix("invG", Nsites, Nsites, invG);

    complex<double>** G = Array2D< complex<double> >(Nsites, Nsites);   
    InvertSymmetricMatrix(Nsites, invG, G);

    //if(n==0) PrintMatrix("G", Nsites, Nsites, G);

    for(int id=0; id<Nsites; id++)
    {  a->r[id].Delta[n] = ii*a->r[id].omega[n] + a->r[id].mu - a->r[id].Sigma[n] - 1.0/G[id][id];
       if (PatchDelta) a->r[id].PatchDelta(AtomicCutoff);
    }

    FreeArray2D< complex<double> >(invG, Nsites);
    FreeArray2D< complex<double> >(G, Nsites);
  }

/*
  if (PrintIntermediate)
  { char FN[300];
    sprintf( FN, "IADHM.U%.3f.T%.3f.it%d", U, T, Iteration);
    r->PrintResult(FN);
  }  
*/  
}
