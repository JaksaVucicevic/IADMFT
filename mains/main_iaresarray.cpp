#include "../source/IAresArray.h"
#include "../source/IAResult.h"
#include "../source/IAGRID.h"
#include <cstdlib>
#include <cstdio>

int main()
{
/*  int N = pow(size,3);
  srand (seed);

  for(int i=0; i<N; i++)
    mus[i] = mu + W*((double) rand()) / ((double) RAND_MAX);
*/

  int N = pow(2,13);
  int Nr = 10;
  
  IAGRID g(N,N, 0.01);

  IAresArray a(Nr, &g);
  a.Set_n(0.5);
  a.Set_mu(1.3);
  for(int i = 0; i<Nr*N; i++) a.totalDelta[i] = 1.0e-5 * i;
  a.ReadTotalDelta();

  printf("a.r[0].n: %.3f a.r[0].mu: %.3f\n",a.r[0].n, a.r[0].mu);

  system("mkdir a");
  a.PrintAll("a/");

  IAresArray b(Nr, &g);
  b.Set_n(0.8);

  a.CopyFrom(b);
  system("mkdir b");
  a.PrintAll("b/");

  b.ReadFromFiles("a/");
  printf("b.r[0].n: %.3f b.r[0].mu: %.3f\n",b.r[0].n, b.r[0].mu);
  system("mkdir a.read");
  b.PrintAll("a.read/");

  return 0;
}

