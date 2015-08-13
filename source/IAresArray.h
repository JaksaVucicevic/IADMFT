#include <complex>
using namespace std;

class IAResult;
class IAGRID;

class IAresArray
{
  public:
    IAresArray(int N, IAGRID* g);
    IAresArray(IAresArray &a); 
    ~IAresArray();

    void Set_n(double n);
    void Set_n0(double n0);
    void Set_mu(double mu);
    void Set_mu0(double mu0);

    double* Get_ns();
    double* Get_mus();
    double* Get_n0s();
    double* Get_mu0s();
    double* Get_Aw0s();
    double* Get_taus();
    double* Get_Zs();
    double* Get_Zs_RASC(complex<double>* Delta=NULL);
    double* Get_renormalized_mus();
    double* Get_ReSigw0s();
    double* Get_pade_Aw0s();
    double* Get_pade_taus();
    double* Get_taus_RASC(complex<double>* Delta=NULL);

/*    double* GetSorted_ns();
    double* GetSorted_mus();
    double* GetSorted_n0s();
    double* GetSorted_mu0s();
    double* GetSorted_Aw0s();
    double* GetSorted_taus();
*/

    IAResult** GetSortedResults();
    double* GetSorted(double* X);

    void Print_ns_and_mus(const char* bareFN);
    void Print_Histogram(const char* FN, double* X, bool Logarithmic=false);
    void Get_Histogram(double* X, int Nbins, double* &x, double* &P, bool Logarithmic=false);

    double Global_n();
    void GetAverageGandSigma(complex<double>* G, complex<double>* Sigma);
    complex<double>* GetGtypical();

    int GetActualN(); //use when reading results from UberMinimal files
    
    IAResult* r;
    complex<double>* totalDelta;
    void WriteTotalDelta();
    void ReadTotalDelta();

    void PrintAll(const char* bareFN);
    void PrintAllShort(const char* bareFN);
    void PrintAllMinimal(const char* bareFN, double cutoff);
    void PrintSortedMinimal(const char* bareFN, double cutoff);
    void PrintAllUberMinimal(const char* bareFN, double cutoff);
    bool ReadFromFiles(const char* bareFN);
    bool ReadFromMinimalFiles(const char* bareFN);
    bool ReadFromUberMinimalFiles(const char* bareFN);
    void CopyFrom(IAresArray &a);
 
    int get_N() {return N;};
  private:
    int N;
};


