# Read application specific variables

include make.inc

# GENERAL PURPOSE MAKEFILE

SRC  := $(wildcard *.f)

MODS := $(wildcard *_mod.f)

OBJS := $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))

OBJMOD := $(MODS:.f=.o)

LIBS :=

COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.f)))
COMMON_OBJS = $(COMMON_SRC:.f=.o)

prefix = .

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O1
  DEBUG = -g -et -Rp
  #DEBUG = -g -et --trace=ALL -Rb
  #DEBUG = -g
  PROFILE = -P
  STATIC = -s
  MODPATH = -p
  FFLAGS = -w

  LIBS = -llapack_f90 -lblas_f90

  MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH)$(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
  DEBUG = -g --chkglobal
  PROFILE =
  STATIC = 
  MODPATH = -M
  FFLAGS =

  LIBS = -llapack_lf95 -lblas_lf95

  MODDIRS = $(MODPATH) . $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = O2 -mp -axW
  DEBUG = -g -check -traceback
  PROFILE = -p
  STATIC =
  MODPATH =
  FFLAGS =

  LIBS = -llapack_intel -lblas_intel

  MODDIRS = $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

#Setup compiler options

OPT = O

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

target: common $(OBJMOD) $(OBJS) $(COMMON_SRC)

common: ;
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -e -C $$subdir ; done 

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS) : $(MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(MODDIRS) $(FFLAGS) $<
