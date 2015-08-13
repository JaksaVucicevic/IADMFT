#include <iostream>
#include "Loop.h"
#include "FFT.h"
using namespace std;

class IASIAM;

class IADHM: public Loop
{
  protected:
    void Defaults();
    string ParamsFN;

    double U;
    double T;
    double W;

    IASIAM* iasiams;

/*   
    double t;
    bool UseBethe; //Delta=t^2 G
                   //... or  ..
    int Nw;        //fill in the non-interacting DOS. and then G_latt= int deps rho_0(eps) (iw+mu-eps-Sigma)^{-1}
    double* w;     //eps
    double* NIDOS; //rho_0
*/

    virtual bool SolveSIAM();
    virtual void CalcDelta();  
   
    virtual void ReleaseMemory();

  public:
    IADHM();
    IADHM(const char* ParamsFN);   
    ~IADHM();

    FFT* fft;

    double** H0;
  
    void SetParams( double U, double W, double T );
    void RandomizeMus( IAresArray &a, double mu, unsigned int seed=0 );
    void GetMus( IAresArray &a, double* mus );
    void SetMus( IAresArray &a, double* mus );
    //void SetNIDOS(int Nw, double* NIDOS, double* w);

    bool PHSymmetricCase;

    bool PatchTailWithAtomicLimit;
    double AtomicCutoff;
    bool PatchDelta;
    
    //void SetUseBethe(bool UseBethe);

    double get_U() { return U; };
    double get_T() { return T; };
};


