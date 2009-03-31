# This is the makefile for LSODE on Linux computers.

FFLAGS = -B100 -O
FC = f90

.f.o: 
	$(FC) $(FFLAGS) -c $*.f

OBJS = \
	daxpy.o \
	dcfode.o \
	dcopy.o \
	ddot.o \
	dewset.o \
	dgbfa.o \
	dgbsl.o \
	dgefa.o \
	dgesl.o \
	dintdy.o \
	dprepj.o \
	dscal.o \
	dsolsy.o \
	dsrcom.o \
	dstode.o \
	dumach.o \
	dvnorm.o \
	idamax.o \
	iumach.o \
	ixsav.o \
	lsode.o \
	xerrwd.o \
	xsetf.o \
	xsetun.o

target: lib

lib: $(OBJS)
	ar -r liblsode.a *.o

clean:
	rm *.o *.a 2> /dev/null

distclean: clean
