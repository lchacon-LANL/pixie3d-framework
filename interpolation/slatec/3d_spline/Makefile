SRC  := $(wildcard *.f)

MODS := $(wildcard *_mod.f)

OBJS := $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))

OBJMOD := $(MODS:.f=.o)

LIBS :=
 
#Define compiler flags

FC = f90
#FC = lf95

ifeq ($(FC),f90)

# Flags for Absoft f90
OPTIMIZATION = -O2 -B101 -cpu:host -s
DEBUG = -g
PROFILE = -P

else

# Flags for Lahey lf95
OPTIMIZATION = -O --tpp
DEBUG = -g --trace
PROFILE =

endif

#Setup compiler options

OPT = O

FFLAGS =

ifneq (,$(findstring O,$(OPT)))
FFLAGS += $(OPTIMIZATION)
endif
ifneq (,$(findstring g,$(OPT)))
FFLAGS += $(DEBUG)
endif
ifneq (,$(findstring p,$(OPT)))
FFLAGS += $(PROFILE)
endif

#Define targets

target: $(OBJMOD) $(OBJS)

clean:
	rm -f *.o *.mod

#Define dependencies

$(OBJS) : $(OBJMOD)

#Define patterns

%.o : %.f
	$(FC) -c $(FFLAGS) $<
