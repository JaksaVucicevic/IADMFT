#include <complex>

using namespace std;

void pade( int M, double* iw, complex<double>* Giw, 
           int N, double* w, complex<double>* Gw );

void PadeToFile( int Mmax, complex<double>* Giw,  double* iw, const char* outputFN, int N, double wmax );

int PadeFromFile(const char* FN, int Mmax, double wmax, int N);

