PROG =	herm

SRCS =	herm.f90 hermite_gauss_module_old.f90

OBJS =	herm.o hermite_gauss_module_old.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O -g -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

herm.o: hermite_gauss_module_old.o
