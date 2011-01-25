# Buildtype configuration

BUILDTYPE?=$(shell echo "`uname -m`-`uname -s`")

ifeq "$(BUILDTYPE)" "i386-Darwin"
  FC=gfortran
#  ifdef SAMR
    HDF5_HOME =${HDF_HOME}
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  =-I$(HDF5_HOME)/include
    PETSC_DIR = /opt/petsc-3.0.0-p9-gcc-4.4.1-debug
    PETSC_ARCH = darwin9.8.0-c-debug
#  endif
endif

# Machine configuration
ifndef HOSTNAME
    HOSTNAME = $(shell hostname)
endif
HOST?=$(HOSTNAME)

MACHINE ?=HaveNoIdea

ifeq ($(HOST),lap0086227)
#    FC=gfortran
    FC=vtf77
    OPT=g
#    HDF5_HOME =/usr
#    HDF5_HOME =$(HDF5)
#    HDF5_HOME =/packages/hdf5-1.8.3/hdf5
    HDF5_HOME =/usr
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  =-I$(HDF5_HOME)/include
    # PETSC_DIR = ${PETSC_DIR}
    # PETSC_ARCH = ${PETSC_ARCH}
    LIBS += -L${PETSC_DIR}/${PETSC_ARCH}/lib -lflapack -lfblas
endif

ifeq ($(HOST),quevedo.ornl.gov)
  MACHINE=quevedo

  FC = gfortran

  ifdef SAMR

    PETSC_DIR        =/opt/petsc-dev
    PETSC_ARCH       =openmpi-gcc-4.4.5-debug
    SAMRAI           =/opt/samrai/samrai-v2.4.4-petsc-dev-openmpi-gcc-4.4.5-debug
    AMRUTILITIES_HOME=/opt/samrutils
    SAMRSOLVERS_HOME =/opt/samrsolvers

  else
    ifdef BOPT
      PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

      FC=gfortran
      ifeq ($(BOPT),O)
        PETSC_ARCH = linux64_openmpi
      else
        PETSC_ARCH = linux64_openmpi_g
      endif

      HDF5 = t

      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include

      ADIOS = t

      ADIOS_HOME = $(HOME)/lib/adios-1.2.1
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml

    else
      LIBS += -L/usr/lib64 -llapack -lblas
      HDF5 = f
      ifeq ($(FC),f90)
        HDF5 = t

        HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
        HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  = $(HDF5_HOME)/lib
        HDF5_INC  =-I$(HDF5_HOME)/include

      endif
      ifeq ($(FC),gfortran)
        HDF5 = t

//         HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/serial
//         HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
//         HDF5_MOD  = $(HDF5_HOME)/lib
//         HDF5_INC  =-I$(HDF5_HOME)/include

        HDF5_HOME =/usr/lib64
        HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  =$(HDF5_HOME)/gfortran/modules
        HDF5_INC  =-I/usr/include

        ADIOS = f
      endif
    endif

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

    VMEC = t
#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

endif

ifeq ($(HOST),gongora.ornl.gov)
  MACHINE=gongora

  FC=gfortran

  LIBS += -L/usr/lib -llapack -lblas

  PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7
  PETSC_ARCH = linux_openmpi

  HDF5 = t
  HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
  HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
  HDF5_MOD  = $(HDF5_HOME)/lib
  HDF5_INC  =-I$(HDF5_HOME)/include

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

  VMEC = t
#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),siesta)
  MACHINE=siesta

  FC=gfortran

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    FC=gfortran
    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi_g
    endif

    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  =-I$(HDF5_HOME)/include
  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f
    ifeq ($(FC),f90)
      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include
      HDF5 = t
    endif
    ifeq ($(FC),gfortran)
      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include
      HDF5 = t
    endif
  endif

  VMEC = t

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),ra22)
  MACHINE=ra22

  FC     = pgf95

  PETSC_ARCH = linux64_openmpi
  PETSC_DIR =/ricercatori/ft/petsc-3.0.0-p7

  override HDF5 = f

  CPPFLAGS += $(PREPROC)RFX
  LIBS     += -llapack -lblas
endif

ifeq ($(HOST),franklin)
  MACHINE=franklin

  FC = pgf95
  HDF5_HOME = $(HDF5_PAR_DIR)
#  HDF5_LIBS:= $(HDF5)
  override HDF5 = t
  BOPT = t
endif

ifeq ($(findstring smoky,$(HOST)),smoky)
  MACHINE=smoky

  FC = pgf95
#  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5 = t
  BOPT = t
endif

ifeq ($(findstring jaguar,$(HOST)),jaguar)
  FC = pgf95

  BOPT = t
  PETSC_DIR ?= $(CRAY_PETSC_DIR)
  PETSC_ARCH ?= cray-xt

  HDF5 = t
  HDF5_HOME = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib

  ADIOS = t
  ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) $(ADIOSREAD_FLIB)

  MACHINE = jaguar
endif

ifeq ($(findstring chester,$(HOST)),chester)
  # SAME AS JAGUAR
  FC = pgf95

  BOPT = t
  PETSC_DIR ?= $(CRAY_PETSC_DIR)
  PETSC_ARCH ?= cray-xt

  HDF5 = t
  HDF5_HOME = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib

  ADIOS = t
  ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) $(ADIOSREAD_FLIB)

  MACHINE = chester
endif

ifeq ($(findstring ewok,$(HOST)),ewok)
  FC = pgf95
  HDF5_DIR ?= /sw/ewok/hdf5/1.6.7/sl5.2_pgi7.2.4_ompi1.2.8
  HDF5_FLIB ?= -module . -module /sw/ewok/hdf5/1.6.7/sl5.2_pgi7.2.4_ompi1.2.8/lib -I/sw/ewok/hdf5/1.6.7/sl5.2_pgi7.2.4_ompi1.2.8/include -L/sw/ewok/hdf5/1.6.7/sl5.2_pgi7.2.4_ompi1.2.8/lib -lhdf5_fortran -lhdf5 -L/sw/ewok/szip/2.1/sl5.2_pgi7.2.4/lib -lsz -lz
  #HDF5_DIR ?= /sw/ewok/hdf5/1.8.2/sl5.2_pgi7.2.4_par
  #HDF5_FLIB ?= -module . -module /sw/ewok/hdf5/1.8.2/sl5.2_pgi7.2.4_par/lib -I/sw/ewok/hdf5/1.8.2/sl5.2_pgi7.2.4_par/include -L/sw/ewok/hdf5/1.8.2/sl5.2_pgi7.2.4_par/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L/sw/ewok/szip/2.1/sl5.2_pgi7.2.4/lib -lsz -lz -lm
  HDF5_HOME = $(HDF5_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib
  HDF5 = t

  #ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_HOME = /ccs/proj/e2e/pnorbert/ADIOS/ewok.noxml
  #ADIOS_LIBS = $(ADIOS_LIB) 
  ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
  ADIOS = t

  OPTIMIZATION = -O2 -Mbounds
  DEBUG = -g -Mbounds -Mchkptr
  BOPT = t
  MACHINE = ewok
endif

ifeq ($(HOST),chenshp.ornl.gov)
    MACHINE=chenshp

    LIBS += -L/usr/lib64/atlas -llapack -lblas
    HDF5=f
endif
