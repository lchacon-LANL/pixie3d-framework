SRC  := $(wildcard *.f)

OBJS := $(SRC:.f=.o)

MODS := $(wildcard *_mod.f)

OBJMOD := $(MODS:.f=.o)

LIBS = 

#FC = lf95
#FFLAGS = -O --tpp
#DEBUG = -g --trace


OPTIMIZATION = -O2 -B101 -cpu:host -s
DEBUG = -g
PROFILE = -P

FC = f90
FFLAGS =  -c $(OPTIMIZATION) $(DEBUG) $(PROFILE)
LDFLAGS = 

target: $(OBJMOD) $(OBJS)

clean:
	rm -f *.o *.mod

$(OBJS) : $(OBJMOD)
