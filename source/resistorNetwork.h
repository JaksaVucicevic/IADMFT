double Resistivity(int Direction, //0: x, 1: y, 2: z, 3: average
                   long long int Nx, 
                   long long int Ny,
                   long long int Nz,
                      long long int BoundaryX,
                      long long int BoundaryY,
                      long long int BoundaryZ,
                   double* rhos,
                   bool SelfNeighbor = false);
