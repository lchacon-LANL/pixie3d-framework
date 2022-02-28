#Define compiler flags

#Defaults

FC  ?= gfortran
OPT ?= O
PREPROC ?= -D

# Flags for Absoft f90
ifeq ($(findstring f90,$(FC)),f90)
  OPTIMIZATION ?= -O3 -cpu:host
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
  LABEL        = absoft
endif

# Flags for Lahey lf95
ifeq ($(findstring lf95,$(FC)),lf95)
  OPTIMIZATION ?= --nap --nchk --fast --npca --nsav --ntrace
  DEBUG        = -g --chk aseu --warn --f95 --trace --trap
#  DEBUG        = -g --chkglobal --warn --f95 --trace
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  VERBOSE      = --verbose
  CPPFLAGS    += $(PREPROC)lf95 
#  FFLAGS      += --quiet --li -Wa,--32
  FFLAGS      += --quiet --li
  LIBS         = -llapackmt -lblasmt
#  LDFLAGS      = -L/usr/lib32 -m32

  ifdef BOPT
    FFLAGS    += --ml cdecl
  endif
  LABEL        = lahey
endif

# Flags for Intel ifort
ifeq ($(findstring ifort,$(FC)),ifort)
  OPTIMIZATION ?= -O3 -g -ip -qopenmp-simd  -align all -march=$(M_ARCH) -fma
#  DEBUG        = -g -check bounds -traceback
  DEBUG        = -g -check noarg_temp_created,bounds,format,pointers,stack -traceback -fpe0 #-warn interfaces
#For threading debugging  DEBUG        = -g -check none -shared-intel #-warn interfaces
  PROFILE      = -O2 -g -p
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)ifort
  FOPENMP      = -qopenmp 
  LABEL        = intel
  LIBS         = -llapack -lblas
endif

# Flags for g95
ifeq ($(findstring g95,$(FC)),g95)
  OPTIMIZATION ?= -O2
  DEBUG        = -g -fbounds-check -ftrace=full
  PROFILE      = -O2 -g -pg
  STATIC       = -fstatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)g95
  FFLAGS      += $(STATIC)
  LABEL        = g95
endif

# Flags for gfortran
ifeq ($(findstring gfortran,$(FC)),gfortran)
  ifeq ($(ARCH),Darwin)
    OPTIMIZATION ?= -Wa,-q -Wa,-Wno-unused-command-line-argument -O3 -g -march=$(M_ARCH) -fopenmp-simd #-fopt-info -ffast-math 
  else
    OPTIMIZATION ?= -O3 -g -march=$(M_ARCH) -fopenmp-simd -no-pie #-fopt-info  -ffast-math 
  endif
  DEBUG        = -g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero -ftrapv #-finit-real=snan #-finit-local-zero #-ggdb3
#  DEBUG        = -g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,denormal -ftrapv #-finit-real=snan #-finit-local-zero #-ggdb3
  PROFILE      = -O2 -pg -g
  STATIC       = -fno-automatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)gfortran
  FFLAGS      += -std=legacy #$(STATIC)
  FOPENMP      = -fopenmp
  LDFLAGS     += -lgfortran 
  LABEL        = gnu
  LIBS         = -llapack -lblas
endif

ifeq ($(findstring vtf77,$(FC)),vtf77)
  OPTIMIZATION ?= -O3 -march=native
  DEBUG        = -ggdb -fbounds-check -fbacktrace -ffpe-trap=invalid,zero -ftrapv
  PROFILE      = -O2 -pg -g
  STATIC       = -fno-automatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)gfortran
  FFLAGS      += #$(STATIC)
  LDFLAGS     += -lgfortran
  LABEL        = gnu
endif

# Flags for Portland Group f90
ifeq ($(findstring pgf95,$(FC)),pgf95)
  OPTIMIZATION ?= -fastsse
#  OPTIMIZATION ?= -fastsse -Mipa=fast
#  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=fp
  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=divz -traceback
  PROFILE      = -pg
  STATIC       = 
  MODFLAG      = -module 
  ADDMODFLAG   = -module 
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)pgf90
  FOPENMP      = -mp
  LABEL        = pgf90
endif

#Assemble compiler options

MODPATH = $(MODFLAG).

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
ifneq (,$(findstring om,$(OPT)))
FFLAGS += $(FOPENMP)
endif

export FFLAGS MODFLAG ADDMODFLAG LDFLAGS PREPROC OPT LABEL CPPFLAGS
