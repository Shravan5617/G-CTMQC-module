FC     = ifort
CFLAG  = -Wall -pedantic -O0 -fopenmp -Wmaybe-uninitialized -lm  ##-fcheck=all -g -fbacktrace -ffpe-trap=invalid
LFLAG  = -llapack -L/usr/local/opt/lapack/lib
LDFLAG = -L. -lpot  

.SUFFIXES : .o .f90 .f
.f90.o:
	$(FC) $(CFLAG) -c $<
.f.o:
	$(FC) $(CFLAG) -c $<

%.o: %.mod

MAIN = main.o

MODS = \
shopping.o \
atomic_masses.o \
classical_evolution.o \
coefficients_evolution.o \
analytical_potentials.o \
electronic_problem.o \
coherence_corrections.o \
kinds.o \
output.o \
trajectories_selection.o \
time_evolution.o \
tools.o \
variables.o \
wigner_distribution.o

all:  main.x

main.x: $(MAIN) $(MODS) $(SUBS)
	$(FC) $(CFLAG) -o $@ $(MAIN) $(MODS) $(SUBS) \
	$(JDLIB) $(WSMPLIB) $(LFLAG) $(LDFLAG)

.PHONY: clean
clean:
	rm -f *.o *.mod *.x

include make.depends
