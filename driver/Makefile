SRC  := $(wildcard *.f)

#MODS := $(wildcard *_mod.f)

MODS := driver_mod.f graphics_mod.f

OBJS := $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))

OBJMOD := $(MODS:.f=.o)

LIBS :=

COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.f)))
COMMON_OBJS = $(COMMON_SRC:.f=.o)

prefix = .

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

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/grid \
	  $(COMMONDIR)/newtonGMRES

#Define targets

target: common $(OBJMOD) $(OBJS) $(COMMON_SRC)

common: ;
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -e -C $$subdir ; \
		cp $$subdir/*.mod $(prefix) >& /dev/null ; done

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS) : $(MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(FFLAGS) $<
