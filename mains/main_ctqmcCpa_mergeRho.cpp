#include <cstdlib>
#include <cstdio>


//=================================== ALL vs T ============================================//


int main(int argc, char* argv [])
{
  if(argc<2) exit(0);
  double W = atof(argv[1]);
//  double T = atof(argv[1]);
//  double U = atof(argv[2]);
//  double mu=U/2.0; 


  // prepare files
  char rhoRASCWFN[300];
  sprintf(rhoRASCWFN,"rho_RASC.W%.3f.new",W);
  FILE* rhoRASCWFile = fopen(rhoRASCWFN,"w");
    

  // ===================== ITERATE ======================== //
  for(double U = 1.9; U<W-1.0; U+=0.5)
  {
      char rhoRASCWUFN[300];
      sprintf(rhoRASCWUFN,"rho_RASC.W%.3f.U%.3f",W,U);
      FILE* rhoRASCWUFile = fopen(rhoRASCWUFN,"r");

      char ln[1000];
      
      while ( fgets (ln , 999 , rhoRASCWUFile) != NULL )
        fprintf(rhoRASCWFile, "%.15le %s", U, ln);    

      fclose(rhoRASCWUFile);

      fprintf(rhoRASCWFile, "\n");

  }
  for(double U = W-1.0; U<W+1.0; U+=0.05)
  {
      char rhoRASCWUFN[300];
      sprintf(rhoRASCWUFN,"rho_RASC.W%.3f.U%.3f",W,U);
      FILE* rhoRASCWUFile = fopen(rhoRASCWUFN,"r");

      char ln[1000];
      
      while ( fgets (ln , 999 , rhoRASCWUFile) != NULL )
        fprintf(rhoRASCWFile, "%.15le %s", U, ln);    

      fclose(rhoRASCWUFile);

      fprintf(rhoRASCWFile, "\n");

  }
 
  fclose(rhoRASCWFile);

  return 0;
}





