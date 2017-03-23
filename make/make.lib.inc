
###################################################
######## Generic package configuration ############
###################################################

# LSODE setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/lsode -llsode

# SLATEC setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/slatec/lib -lslatec
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/slatec/lib

# FPA setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/fpa/lib -lfpa
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/fpa/lib

# SDC setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/sdc/lib -lsdc
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/sdc/lib

# ADIOS setup

ifeq ($(ADIOS),t)
  CONTRIBLIBS += $(ADIOS_LIBS)
  CPPFLAGS    += $(PREPROC)adios $(PREPROC)ADIOS_VER=$(ADIOS_VER)
  ifdef ADIOS_BUFFER_MB
    CPPFLAGS  += $(PREPROC)ADIOS_BUFFER_MB=$(ADIOS_BUFFER_MB)
  endif
  MODPATH     += $(ADDMODFLAG)$(ADIOS_INC)
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS += $(HDF5_LIBS) 
  CPPFLAGS    += $(PREPROC)hdf5 $(PREPROC)H5_USE_16_API -I$(HDF5_INC)
  MODPATH     += $(ADDMODFLAG)$(HDF5_MOD)
endif

# ARPACK setup

ifeq ($(ARPACK),t)
  CONTRIBLIBS += $(ARPACK_LIBS)
  CPPFLAGS    += $(PREPROC)arpack
endif

# PETSC setup

ifdef BOPT
  include ${PETSC_DIR}/conf/base

  CPPFLAGS += $(PREPROC)petsc -I$(PETSC_DIR)/include -I${PETSC_DIR}/$(PETSC_ARCH)/include

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif
endif

#SAMRAI setup

ifdef SAMR
  ifdef SAMRSOLVERS_HOME
     include ${SAMRSOLVERS_HOME}/Makefile.config
  else
     $(error error SAMRSOLVERS_HOME must be set)
  endif

  OBJECT=${SAMRAI}

  CPPFLAGS += $(PREPROC)samrai
endif

# RNG setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/rng/lib -lrng

# BTRIDIAG (block tridiagonal solver) setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/btridiag/lib -lbtrid
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/btridiag/lib

# PTRIDIAG (parallel tridiagonal solver) setup

ifeq ($(PTRID),t)
  CONTRIBLIBS += -L$(COMMONDIR)/contrib/ptridiag/lib -lptrid
endif

# SPRNG setup

ifeq ($(SPRNG),t)
  CONTRIBLIBS += -L$(SPRNG_LIBS) 
  CPPFLAGS    += $(PREPROC)_SPRNG -I$(SPRNG_INC) -DHAVE_CONFIG_H -DPOINTER_SIZE=8 -DLONG64=long -DINTEGER_STAR_8
#  MODPATH     += $(ADDMODFLAG)$(SPRNG_MOD)
endif
