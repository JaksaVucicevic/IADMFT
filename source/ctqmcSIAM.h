#include <string>

class IAResult;

using namespace std;

class ctqmcSIAM
{
  public:
    ctqmcSIAM();
    //ctqmcSIAM(const ctqmcSIAM &s);
    ~ctqmcSIAM();

    void Defaults();

    double GetOccupancy();

    void PrintParamsFile();
    void PrintCixFile();
    void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X);
    void ReadXY(const char* FileName, int N, complex<double>* Y, double* X);


    bool Run(IAResult* r);

    
    //------------ options---------------//

    bool DoPrintParamsFile;
    bool UseSmartNom;
    double FreqCutoff;

    bool PatchTailWithAtomicLimit;
    double AtomicCutoff;
    
    //------------ params --------------//

    string execName;
    string runFolderName; //full path!
    bool runParallel;

    string inputf;      // filename of the parameter file
    string fDelta; // input bath function Delta(iom)
    string fcix;   // input file with atom eigenvalues and eigenvectors (exact diagonalization results)
    string fGf;    // output file for green's function G(iom)
    string fSig;   // output file for self-energy Sig(iom)
    double mu;           // chemical potential
    double beta;         // inverse temperature
    double U;           // Coulomb repulsion (should be used only in single site DMFT)
    unsigned long int M;     // number of all qmc steps
    int    Ntau;         // number of imaginary time points for Delta(tau)
    int    Nmax;        // maximum number of kinks
    int    nom ;         // number of frequency points sample in qmc
    int    nomb ;          // number of bosonic frequencies for susceptibilities
    double PChangeOrder;     // probability to change the order (as compared to move a kink)
    int    tsample;          // how often to record measurements (each tsample steps)
    int    warmup;      // how many qmc steps should be ignored
    int CleanUpdate;      // clean update is sometimes useful
    double minM;       // trace shuld always be larger than this minimal value when accepting the step
    double minD;       // determinant of hybrodization should always be larger than this number when accepting the step
    int    Ncout;       // screen output for basic information every Ncout steps
    int    Naver;     // much more information about simulation every Naver steps
    double TwoKinks;          // probability to add two kinks
    int    GlobalFlip;         // global flip is tried after GlobalFlip qmc steps
    double treshold;       // energy to throw away atomic states when creating new.cix for next iteration
    int  SampleGtau;          // How often to sample for Gtau (if at all)
    int SampleVertex;          // How often to sample two particle vertex
    int  ReduceComm;           // only part of Gtau is transfered between nodes (if slow communication)
    int    Ncorrect;          // baths with index higher than Ncorrect should not be added a high-frequency tail (if -1, all tails are added)
    int         aom;           // number of frequency points to find Sigma(om_last_sampled)
    int         som;           // number of frequency points to find Susc(om_last_sampled)
    int    PreciseP;           // computes probabilities more precisely
    double   sderiv;         // maximum mismatch when matching high and low frequency of imaginary part of Sigmas
    double minDeltat;        // Delta(tau) is sometimes not causal due to numerical error. In this case we set it to small value minDeltat.
    bool SampleSusc;       // If spin and charge dynamic susceptibility should be sampled during simulation
    int        nomv;       // number of fermionic frequencies for computing vertex
    int         nOm;           // number of bosonic frequencies for vertex calculation

};
