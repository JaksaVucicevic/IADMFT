#include <complex>
using namespace std;

class IAGRID;
class IBZ;

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
    complex<double>* Delta;	//Bath green's function
    complex<double>* G0;	//Bath+non-interacting impurity Green's function
    complex<double>* SOCSigma;	//Second order in the interaction expansion in sigma Sigma
    complex<double>* Sigma;	//Total Sigma, + corrections to interpolate between exact limiting cases
    complex<double>* G;		//Impurity Green's function

    void PrintResult(const char* FN);
    bool ReadFromFile(const char* FN);
    void CopyFrom(const IAResult &r);

    complex<double> OpticalConductivity(int n, IBZ* ibz);
    complex<double> Lambda(int n, IBZ* ibz);

    complex<double> OpticalConductivity(int n, double (*v)(double)=NULL);
    complex<double> Lambda(int n, double (*v)(double)=NULL);
        
    //double TriangularConductivity(double T, int Nkx, int Nky, int Nnu, const char * integrandFN=NULL);

  private:
    void Initialize(IAGRID* g);
    void ReleaseMemory();
};
