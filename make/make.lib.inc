
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
  CPPFLAGS    += $(PREPROC)adios $(PREPROC)ADIOS_VER=$(ADIOS_VER)  $(PREPROC)ADIOS_FILE=\"record-$(ADIOS_VER).bp\"
  ifdef ADIOS_BUFFER_MB
    CPPFLAGS  += $(PREPROC)ADIOS_BUFFER_MB=$(ADIOS_BUFFER_MB)
  endif
  ifdef ADIOS_METHOD
    CPPFLAGS  += $(PREPROC)ADIOS_METHOD=\"$(ADIOS_METHOD)\"
  endif
  MODPATH     += $(ADDMODFLAG)$(ADIOS_HOME)/include
  export ADIOS_CONFIG = adios1.10.0_config.xml
  export ADIOS_VER
  export ADIOS
endif

# ADIOS2 setup

ifeq ($(ADIOS2),t)
  CONTRIBLIBS += $(ADIOS2_LIBS)
  ## jyc: ADIOS2 should be capital
  CPPFLAGS    += $(PREPROC)ADIOS2 $(PREPROC)ADIOS2_VER=$(ADIOS2_VER)  $(PREPROC)ADIOS2_FILE=\"record-$(ADIOS2_VER).bp\"
  ifdef ADIOS2_METHOD
    CPPFLAGS  += $(PREPROC)ADIOS2_METHOD=\"$(ADIOS2_METHOD)\"
  endif
  MODPATH     += $(ADIOS2_MOD)
  export ADIOS_CONFIG = adios2.5.0_config.xml
  export ADIOS_VER = $(ADIOS2_VER)
  export ADIOS2
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS += $(HDF5_LIBS) 
  CPPFLAGS    += $(PREPROC)hdf5 -I$(HDF5_INC)
  MODPATH     += $(ADDMODFLAG)$(HDF5_MOD)
endif

# ARPACK setup

ifeq ($(ARPACK),t)
  CONTRIBLIBS += $(ARPACK_LIBS)
  CPPFLAGS    += $(PREPROC)arpack

  export ARPACK
endif

# PETSC setup

ifdef BOPT
#  include ${PETSC_DIR}/conf/base
  include ${PETSC_DIR}/lib/petsc/conf/variables

#  To access the PETSc variables for the build, including compilers, compiler flags, libraries etc but
#  manage the build rules yourself (rarely needed) comment out the next lines
#
#  include ${PETSC_DIR}/lib/petsc/conf/rules
#  include ${PETSC_DIR}/lib/petsc/conf/test

  CPPFLAGS += $(PREPROC)petsc -I$(PETSC_DIR)/include -I${PETSC_DIR}/$(PETSC_ARCH)/include

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif

  export PETSC_DIR PETSC_ARCH FLINKER PETSC_SNES_LIB BOPT SNES_OPT
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

  export CPPFLAGS_EXTRA \
         CXXFLAGS_EXTRA LDFLAGS_EXTRA LDLIBS_EXTRA SAMR SAMRAI LIBSAMRAI3D LIBSAMRAI \
         AMRUTILITIES_HOME SAMRSOLVERS_HOME
endif

# RNG setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/rng/lib -lrng
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/rng/lib

# BTRIDIAG (block tridiagonal solver) setup

CONTRIBLIBS += -L$(COMMONDIR)/contrib/btridiag/lib -lbtrid
MODPATH     += $(ADDMODFLAG)$(COMMONDIR)/contrib/btridiag/lib

# PTRIDIAG (parallel tridiagonal solver) setup

ifdef BOPT
  CONTRIBLIBS += -L$(COMMONDIR)/contrib/ptridiag/lib -lptrid
endif

# SPRNG setup

ifeq ($(SPRNG),t)
  CONTRIBLIBS += -L$(SPRNG_LIBS) 
  CPPFLAGS    += $(PREPROC)_SPRNG -I$(SPRNG_INC) -DHAVE_CONFIG_H -DPOINTER_SIZE=8 -DLONG64=long -DINTEGER_STAR_8
#  MODPATH     += $(ADDMODFLAG)$(SPRNG_MOD)
endif
