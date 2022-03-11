# This is the makefile for LSODE on Linux computers.

FFLAGS = -B100 -O
FC = f90

.f.o: 
	$(FC) $(FFLAGS) -c $*.f
.F.o: 
	$(FC) $(FFLAGS) -c $*.F

OBJS := $(patsubst %.f,%.o,$(wildcard *.f)) $(patsubst %.F,%.o,$(wildcard *.F))

target: lib

lib: $(OBJS)
	ar -r liblsode.a *.o

clean:
	-@rm *.o *.a 2> /dev/null

distclean: clean
