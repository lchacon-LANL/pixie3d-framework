# This is the makefile for LSODE on Linux computers.

FFLAGS = -B100 -O
FC = f90

.f.o: 
	$(FC) $(FFLAGS) -c $*.f

OBJS := $(patsubst %.f,%.o,$(wildcard *.f))

target: lib

lib: $(OBJS)
	ar -r liblsode.a *.o

clean: ;
	-rm -f *.o *.a 

distclean: clean
