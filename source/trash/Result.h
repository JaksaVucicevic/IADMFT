#include <complex>
using namespace std;

class IAGRID;

class IAResult
{
  public:
    IAResult(IAGRID* g);
    IAResult(const IAResult &r);
    ~IAResult();

    void Reset();
    void Reset(IAGRID* g);
 
    IAGRID* iagrid;

    double n;
    double n0;
    double mu;
    double mu0;

    double* tau;
    complex<double>* G0_tau;	//auxillary Green's function
    complex<double>* SOCSigma_tau; 

    double* omega;		//omega grid
    complex<double>* SOCSigma;	//Second order contribution in sigma
    complex<double>* Sigma;	//Sigma interpolating between exact limiting cases
    complex<double>* G;		//Greens function on real axis
    complex<double>* Delta;	//Bath
    complex<double>* G0;	//auxillary Green's function

    double* w;
    double* NIDOS;		//non-interacting density of states

    void PrintResult(const char* FN);
    bool ReadFromFile(const char* FN);
    void CopyFrom(const IAResult &r);
    void PrintModel(double U, double T);
        
    //double TriangularConductivity(double T, int Nkx, int Nky, int Nnu, const char * integrandFN=NULL);
    //void ChargeSusceptibility(double T, double &chi1, double &chi3);

  private:
    void Initialize(IAGRID* g);
    void ReleaseMemory();
};
