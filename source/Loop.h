#include <iostream>

class IAResult;
class IAresArray;
class IAGRID;
class LambdaCalculator;

using namespace std;

class Loop
{
  protected:
    void Defaults();
    string ParamsFN;

    int Iteration;

    IAResult* r;
    IAresArray* a;

    IAGRID* iagrid;

    int N;
    int Nsites;

    //---- Broyden options ----//
    bool UseBroyden;
    bool ForceBroyden;
    double BroydenStartDiff;

    //---- Mixer Options ------//
    int NtoMix;
    int * Coefs;

    //---- Loop Options -------//
    int MAX_ITS;
    int MIN_ITS;
    double Accr;

    //---- PrintOut/Debugging optins----//
    bool PrintIntermediate;
    bool HaltOnIterations;
    bool ForceSymmetry;
    
    //---- Functions to be overridden---//
    virtual bool SolveSIAM();
    virtual void CalcDelta();

    virtual void ReleaseMemory();

  public:
    Loop();
    Loop(const char* ParamsFN);
    ~Loop();
    
    bool Run(IAResult* r);
    bool Run(IAresArray* a);
    
    void SetForceSymmetry(bool FS);
    void SetGrid(IAGRID* g);
    void SetMixerOptions(int NtoMix, const int * Coefs);
    void SetBroydenOptions(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff);
    void SetLoopOptions(int MAX_ITS, double Accr);
    void SetPrintOutOptions(bool PrintIntermediate, bool HaltOnIterations);

    //------------------------------//
    LambdaCalculator* LC;
};
