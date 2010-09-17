#Define compiler flags

#Defaults

FC  ?= gfortran
OPT ?= O
PREPROC ?= -D

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O3 -cpu:host
  DEBUG        = -g
#  DEBUG        = -g -en -et -trap=DIVBYZERO,INVALID
  PROFILE      = -P
  STATIC       = -s
  MODFLAG      = -p
  ADDMODFLAG   = -p
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)absoft_ 

  FFLAGS      += -w -YEXT_NAMES=LCS -YEXT_SFX=_ -YCFRL=1 -m64 -lU77
  LDFLAGS      = -lU77 
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = --nap --nchk --fast --npca --nsav --ntrace
  DEBUG        = -g --chk aseu --warn --f95 --trace
#  DEBUG        = -g --chkglobal --warn --f95 --trace
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  VERBOSE      = --verbose
  CPPFLAGS    += $(PREPROC)lf95 
  FFLAGS      += --quiet

  ifdef BOPT
    FFLAGS    += --ml cdecl
  endif
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = -O2 -mp -axW
  DEBUG        = -g -check all -traceback
#  DEBUG        = -g -check -traceback
  PROFILE      = -p
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)ifort
endif

# Flags for g95
ifeq ($(FC),g95)
  OPTIMIZATION = -O2
  DEBUG        = -g -fbounds-check -ftrace=full
  PROFILE      = -pg
  STATIC       = -fstatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)g95
  FFLAGS      += $(STATIC)
endif

# Flags for gfortran
ifeq ($(FC),gfortran)
  OPTIMIZATION = -O3 -march=native
  DEBUG        = -ggdb -fbounds-check -fbacktrace -ffpe-trap=invalid,zero
  PROFILE      = -pg
  STATIC       = -fno-automatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)gfortran
  FFLAGS      += #$(STATIC) 
endif

# Flags for Portland Group f90
ifeq ($(FC),pgf95)
  OPTIMIZATION = -fastsse -Mipa=fast
  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=fp
  PROFILE      = -pg
  STATIC       = 
  MODFLAG      = -module 
  ADDMODFLAG   = -module 
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)pgf90
endif

#Assemble compiler options

ifneq (,$(findstring O,$(OPT)))
FFLAGS += $(OPTIMIZATION)
endif
ifneq (,$(findstring g,$(OPT)))
FFLAGS += $(DEBUG)
endif
ifneq (,$(findstring p,$(OPT)))
FFLAGS += $(PROFILE)
endif
ifneq (,$(findstring s,$(OPT)))
FFLAGS += $(STATIC)
endif
ifneq (,$(findstring v,$(OPT)))
FFLAGS += $(VERBOSE)
endif
