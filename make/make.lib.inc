
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
   include ${SAMRAI}/config/Makefile.config

   PDIM = 3
   OBJECT=${SAMRAI}
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)

   TARGET = samrai
   CPPFLAGS += -Dsamrai

   SUBDIRS += common/driver-samrai

   CPPFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
   CPPFLAGS_EXTRA += -I${SAMRSOLVERS_HOME}/include
   CXXFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
   CXXFLAGS_EXTRA += -I${SAMRSOLVERS_HOME}/include
   LDFLAGS_EXTRA += -L${AMRUTILITIES_HOME}/lib   
   LDLIBS_EXTRA += ${AMRUTILITIES_HOME}/lib/libAMRUtils3d.a
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/liblinearops3d.a
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libmlsolvers3d.a
   LDFLAGS_EXTRA += -L${SAMRSOLVERS_HOME}/lib   
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libtimeintegrators3d.a
   #LDLIBS_EXTRA += -lAMRUtils${PDIM}d
endif
