
###################################################
######## Generic package configuration ############
###################################################

# LSODE setup

CONTRIBLIBS += -L$(PWD)/common/contrib/lsode -llsode

# SLATEC setup

CONTRIBLIBS += -L$(PWD)/common/contrib/slatec/lib -lslatec
MODPATH     += $(ADDMODFLAG)$(PWD)/common/contrib/slatec/lib

# ADIOS setup

ifeq ($(ADIOS),t)
  CONTRIBLIBS += $(ADIOS_LIBS)
  CPPFLAGS    += $(PREPROC)adios $(PREPROC)ADIOS_VER=$(ADIOS_VER) -I$(ADIOS_HOME)/include
#  MODPATH     += $(ADDMODFLAG)$(ADIOS_HOME)/include
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS += $(HDF5_LIBS) 
  CPPFLAGS    += $(PREPROC)hdf5 $(PREPROC)H5_USE_16_API $(HDF5_INC)
  MODPATH     += $(ADDMODFLAG)$(HDF5_MOD)
endif

# ARPACK setup

ifeq ($(ARPACK),t)
  CONTRIBLIBS += $(ARPACK_LIBS)
  CPPFLAGS    += $(PREPROC)arpack
endif

# FPA setup

ifeq ($(FPA),t)
  CONTRIBLIBS += -L$(PWD)/common/contrib/fpa/lib -lfpa
  CPPFLAGS    += $(PREPROC)FPA
  MODPATH     += $(ADDMODFLAG)$(PWD)/common/contrib/fpa/lib
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