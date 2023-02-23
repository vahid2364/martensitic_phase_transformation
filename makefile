FC=ifort
FCFLAGS=-O -g -traceback -qmkl /sw/eb/sw/imkl/2021.4.0/mkl/2021.4.0/include/mkl_dfti.f90
FLFLAGS=

SOURCES=Allen_cahn_2D.f90
OBJECTS=$(SOURCES:.f90=.o)
TARGET=pfm.out

compile: 
	$(FC) $(FCFLAGS) -c $(SOURCES)

run: 
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJECTS) 

clean:	
	@rm -rf *.o *.mod *.out
