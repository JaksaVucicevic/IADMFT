# C++ compiler
Cpp = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc #/opt/intel/bin/icpc

main =main_iachm

#mpiCC
# source path
SP = ../source
#object files path
OP = ../obj
# executable path
RP = ../bin

# -----------  SERIAL ------------------------------------#

mpiCC = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc
FLAGS =  -lgmp -std c++11 -static-intel

# ----------- THREADED -----------------------------------# 

#mpiCC = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc
#FLAGS =  -D_OMP -openmp #-static

# -----------  HYBRID ------------------------------------#

#mpiCC = /opt/openmpi-1.6.3/bin/mpicxx
#FLAGS =  -D_MPI -D_OMP -openmp -static-intel #-fast

#---------------------------------------------------------#

LIBS =# use this if needed 

all : $(OP)/$(main).o  $(OP)/IACHM.o  $(OP)/Loop.o $(OP)/IASIAM.o $(OP)/FFT.o $(OP)/LambdaCalculator.o $(OP)/IAResult.o $(OP)/IAGRID.o $(OP)/Input.o $(OP)/Broyden.o $(SP)/Broyden.h $(OP)/IBZ.o $(SP)/Mixer.h $(OP)/routines.o $(OP)/GRID.o $(OP)/nrutil.o $(OP)/pade.o
	$(mpiCC) $(FLAGS) -o $(RP)/$(main) $(LIBS) $(OP)/$(main).o  $(OP)/IACHM.o $(OP)/Loop.o $(OP)/IASIAM.o $(OP)/FFT.o $(OP)/LambdaCalculator.o $(OP)/IAResult.o $(OP)/IAGRID.o $(OP)/Input.o $(OP)/Broyden.o $(OP)/IBZ.o $(OP)/routines.o $(OP)/GRID.o $(OP)/nrutil.o $(OP)/pade.o

# main program
$(OP)/$(main).o : $(main).cpp $(SP)/IACHM.h $(SP)/IASIAM.h $(SP)/IAResult.h $(SP)/IAGRID.h $(SP)/FFT.h $(SP)/pade.h
	$(mpiCC) $(FLAGS) -c -o $@ $(main).cpp

# CHM
$(OP)/IACHM.o : $(SP)/IACHM.cpp $(SP)/IACHM.h $(SP)/Loop.h $(SP)/IASIAM.h $(SP)/IAResult.h $(SP)/IAGRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IACHM.cpp

# Loop (base class for CHM and TMT)
$(OP)/Loop.o : $(SP)/Loop.h $(SP)/Loop.cpp $(SP)/IAResult.h $(SP)/IAGRID.h $(SP)/Input.h $(SP)/Mixer.h $(SP)/Broyden.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Loop.cpp

# IASIAM
$(OP)/IASIAM.o : $(SP)/IASIAM.cpp $(SP)/IASIAM.h $(SP)/Broyden.h $(SP)/IAResult.h $(SP)/IAGRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IASIAM.cpp

# FFT 
$(OP)/FFT.o : $(SP)/FFT.h $(SP)/FFT.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/FFT.cpp

# LambdaCalculator
$(OP)/LambdaCalculator.o : $(SP)/LambdaCalculator.cpp $(SP)/LambdaCalculator.h $(SP)/Input.h $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/LambdaCalculator.cpp

# IAResult
$(OP)/IAResult.o : $(SP)/IAResult.cpp $(SP)/IAResult.h $(SP)/IAGRID.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IAResult.cpp

# IAGrid utility for initializing omega grids and provides all grid dependent routines
$(OP)/IAGRID.o : $(SP)/IAGRID.cpp $(SP)/IAGRID.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IAGRID.cpp

# Input class used for reading files with parameters
$(OP)/Input.o : $(SP)/Input.cpp $(SP)/Input.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Input.cpp

# Broyden solver for systems of equations
$(OP)/Broyden.o : $(SP)/Broyden.h $(SP)/Broyden.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Broyden.cpp

# IBZ
$(OP)/IBZ.o : $(SP)/IBZ.h $(SP)/IBZ.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IBZ.cpp

# contains some constants and useful numerical routines
$(OP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# IAGrid utility for initializing omega grids and provides all grid dependent routines
$(OP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# numerical routines from NumRec
$(OP)/nrutil.o : $(SP)/nrutil.h $(SP)/nrutil.c
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/nrutil.c

# numerical routines from NumRec
$(OP)/pade.o : $(SP)/pade.h $(SP)/pade.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/pade.cpp

# clean all object and exec files
clean :
	rm -f $(RP)/$(main) $(OP)/*.o
