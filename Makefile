PROG =	herm

SRCS =	herm.f90 hermite_gauss_module_old.f90 hermite_gauss_module.f90

OBJS =	herm.o hermite_gauss_module_old.o hermite_gauss_module.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O -g -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan -fopenmp -ffpe-trap=zero,overflow,invalid
#F90FLAGS = -O -g -std=f2008 -Wall -Wextra -finit-real=snan -fopenmp 
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

herm.o: hermite_gauss_module_old.o hermite_gauss_module.o
