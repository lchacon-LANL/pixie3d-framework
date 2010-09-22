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

HOST?=$(HOSTNAME)
MACHINE ?=HaveNoIdea

ifeq ($(HOST),quevedo.ornl.gov)
  MACHINE=quevedo
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

      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include

//       HDF5_HOME =/usr/lib64
//       HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
//       HDF5_MOD  =/usr/lib64/gfortran/modules
//       HDF5_INC  =-I/usr/include

      ADIOS = f
    endif
  endif

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include

  VMEC = t
#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
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
  PETSC_ARCH = linux64_openmpi
  PETSC_DIR =/ricercatori/ft/petsc-3.0.0-p7
  override HDF5 = f
  CPPFLAGS += $(PREPROC)RFX
  FC        = pgf95
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
#  ifeq ($(PE_ENV),GNU) 
#    FC = gfortran
#  endif
#  ifeq ($(PE_ENV),PGI)
#    FC = pgf95
#  endif

#  HDF5 = t
##  HDF5_HOME = $(HDF5_PAR_DIR)
#  HDF5_LIBS = $(HDF5_FLIB)

#  BOPT = t
#  PETSC_DIR ?= /opt/petsc/2.3.3/real/${PE_ENV}/${XTPE_COMPILE_TARGET}
#  PETSC_ARCH ?= cray-xt

  FC = pgf95
  #FC = ddtmem-ftn
  HDF5_HOME = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
  #HDF5_LIBS = $(HDF5_FLIB)
  HDF5_LIBS = -L/sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi7.2.3/lib -lsz -lz
  #HDF5_HOME = /sw/xt5/hdf5/1.6.8/cnl2.2_pgi9.0.4
  #HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi7.2.3/lib -lsz -lz
  #HDF5_HOME = /opt/cray/hdf5-parallel/1.8.3.1/hdf5-parallel-pgi
  #HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi8.0.3/lib -lsz -lz
  HDF5_MOD  = $(HDF5_HOME)/lib
  HDF5 = t

  #ADIOS_HOME = $(ADIOS_DIR)
  #ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_LIB) $(ADIOSREAD_LIB)
  ADIOS_HOME = /ccs/proj/e2e/pnorbert/ADIOS/jaguar.noxml
  #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladios -ladios_phdf5 -ladiosread -lmxml
  ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
  ADIOS = t

  OPTIMIZATION = -O2 -Mbounds
  DEBUG = -g -Mbounds -Mchkptr
  #DEBUG = -g -Wl,-z,muldefs -L/sw/xt5/ddt/2.5.1-12323/sles10.1_binary/lib/64 -ldmalloc
  BOPT = t
  MACHINE = jaguar
endif

ifeq ($(findstring chester,$(HOST)),chester)
  # SAME AS JAGUAR
  FC = pgf95
  HDF5_HOME = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
  #HDF5_LIBS = $(HDF5_FLIB)
  HDF5_LIBS = -L/sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi7.2.3/lib -lsz -lz
  #HDF5_HOME = /sw/xt5/hdf5/1.6.8/cnl2.2_pgi9.0.4
  #HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi7.2.3/lib -lsz -lz
  #HDF5_HOME = /opt/cray/hdf5-parallel/1.8.3.1/hdf5-parallel-pgi
  #HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -L/sw/xt5/szip/2.1/sles10.1_pgi8.0.3/lib -lsz -lz
  HDF5_MOD  = $(HDF5_HOME)/lib
  HDF5 = t

  #ADIOS_HOME = $(ADIOS_DIR)
  #ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) $(ADIOSREAD_FLIB)
  ADIOS_HOME = /ccs/proj/e2e/pnorbert/ADIOS/jaguar.noxml
  #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladios -ladios_phdf5 -ladiosread -lmxml
  ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
  ADIOS = t

  OPTIMIZATION = -O2 -Mbounds
  DEBUG = -g -Mbounds -Mchkptr
  BOPT = t
  MACHINE = jaguar
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
