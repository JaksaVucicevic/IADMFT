#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "mkl_invertMatrix.h"

#include <mkl.h>
#include <mkl_types.h>

using namespace std;

void Invert_Matrix_ge(complex<double>* invA, complex<double>* A, int N)
{
  int LWORK=4*N*N;
  int *permutations;
  MKL_Complex16 *WORK, *tempA;
  tempA = new MKL_Complex16 [N*N];
  permutations = new int[2*N];
  WORK = new MKL_Complex16 [4*N*N];
  int *IPIV;
  IPIV = new int[N];
  int INFO=0;

  for (int l=0; l<N*N; ++l)
  {  tempA[l].real=real(A[l]);
     tempA[l].imag=imag(A[l]);
  }
  zgetrf_( &N, &N, tempA , &N, IPIV, &INFO );
  if (INFO != 0)
  {
    cout << "ComplexMatrixInverse: Error at zhetrf INFO = " << INFO; exit(0);
  }
  zgetri_( &N, tempA , &N, IPIV, WORK, &LWORK, &INFO );
  for (int l=0; l<N*N; ++l)
    invA[l]=complex<double>(tempA[l].real,tempA[l].imag);
  if (INFO != 0)
  {
    cout << "ComplexMatrixInverse: Error at zsytri  \n"; exit(0);
  }

  delete [] WORK;
  delete [] tempA;
  delete [] permutations;
  delete [] IPIV;

}

void InvertSymmetricMatrix(int N, complex<double>** A,  complex<double>** invA)
{
  complex<double>* X = new complex<double>[N*N];
  complex<double>* invX = new complex<double>[N*N];
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    X[i*N+j] = A[i][j];

  Invert_Matrix_ge(invX, X, N);

  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    invA[i][j] = invX[i*N+j];

  delete [] X;
  delete [] invX;
}

/*
void PrintMatrix(const char* FileName, int N, int M, complex<double>** A)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
  { for (int j=0; j<M; j++)
      fprintf(f,"%.3f\t", real(A[i][j]));
    fprintf(f,"\n");
  }
  fprintf(f,"\n\n");
  for (int i=0; i<N; i++)
  { for (int j=0; j<M; j++)
      fprintf(f,"%.3f\t", imag(A[i][j]));
    fprintf(f,"\n");
  }

  fclose(f);
}


// TEST-------------------------------------
int main()
{
  int N=5;
  complex<double>** M = new complex<double>*[N];
  complex<double>** invM = new complex<double>*[N];
  for(int i=0; i<N; i++)
  { M[i] = new complex<double>[N];
    invM[i] = new complex<double>[N];  
  }

  for(int i=0; i<N; i++)
  for(int j=i; j<N; j++)
  {  M[i][j] = 1+i+j + (i-j)*(i+j);
     M[j][i] = M[i][j];
  }
  

  PrintMatrix("M", N, N, M);

  complex<double>* A = new complex<double>[N*N];
  complex<double>* invA = new complex<double>[N*N];
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    A[i*N+j] = M[i][j];

  Invert_Matrix_ge(invA, A, N);

  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    invM[i][j] = invA[i*N+j];

  PrintMatrix("invM", N, N, invM);

  return 0;
}
*/
