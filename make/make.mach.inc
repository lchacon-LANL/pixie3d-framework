# Buildtype configuration

BUILDTYPE?=$(shell echo "`uname -m`-`uname -s`")

ifeq "$(BUILDTYPE)" "i386-Darwin"
  FC=gfortran
#  ifdef SAMR
    HDF5_HOME =${HDF_HOME}
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  = $(HDF5_HOME)/include
    PETSC_DIR = /opt/petsc-3.0.0-p9-gcc-4.4.1-debug
    PETSC_ARCH = darwin9.8.0-c-debug
#  endif
endif

# Machine configuration

HOST ?= $(shell hostname -s)

ARCH ?= $(shell uname -s)

CPPFLAGS += $(PREPROC)$(ARCH)

############
# DEFAULTS #
############

FC=gfortran
HDF5 = f
ADIOS= f

###############################
# lap0086227 (Mark Berrill's) #
###############################

ifeq ($(HOST),lap0086227)
    FC=gfortran
    OPT=g
    HDF5 = t
    ifndef SAMR
        HDF5_HOME =/usr
        LIBS += -L/packages/lapack-3.2.2/ -llapack -lblas
    endif
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  =$(HDF5_HOME)/include
    LIBS += -L/packages/acml-4-4-0/gfortran64/lib -llapack -lblas -lgfortran
endif

###############################
# cray2 (CSU cray cluster)    #
###############################

ifeq ($(HOST),cray2)
    FC=gfortran
    OPT=g
    HDF5 = f
    SAMRSOLVERS_HOME=/home/mberrill/lustrefs/samrsolvers-v3/build/optimized
endif

###############################
# vayu (Bobby's desktop)      #
###############################

ifeq ($(HOST),vayu)
    FC=gfortran
    #FC=mpif90
    #CC=mpicc
    #CXX=mpic++
    OPT=O
    HDF5 = t
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  = $(HDF5_HOME)/include
    LIBS += -L${ACML_DIR} -llapack -lblas
endif

################
# jetta (LANL) #
################

ifeq ($(HOST),jetta)

  FC = gfortran

  ifdef SAMR

    CC = mpicc
    CXX = mpicxx
    F77 = gfortran
    override PETSC_DIR  =
#    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/debug
    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/optimized
    HDF5 = f
    HDF5_LIBS = 
    LDFLAGS    += -lgfortran

  else
    ifdef BOPT
      PETSC_DIR ?=$(HOME)/lib/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
        PETSC_ARCH = $(PETSC_ARCH_O)
      else
        PETSC_ARCH = $(PETSC_ARCH_g)
      endif

      HDF5 = t

      #Determined via modules
      HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(PHDF5_LIB)
      HDF5_INC  = $(PHDF5_INC)

      ADIOS = t

      #Determined via modules
      ADIOS_LIBS = -I$(ADIOS_INC) -L$(ADIOS_LIB) -ladiosf -ladiosreadf -lmxml

    else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = t

      #Determined via modules
      HDF5_LIBS =-L$(SHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(SHDF5_LIB)
      HDF5_INC  = $(SHDF5_INC)

      ADIOS = f
    endif

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

endif

################
# mavis (LANL) #
################

ifeq ($(HOST),mavis)

  FC = gfortran

  ifdef SAMR

    CC = mpicc
    CXX = mpicxx
    F77 = gfortran
    override PETSC_DIR  =
#    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/debug
    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/optimized
    HDF5 = f
    HDF5_LIBS =
    LDFLAGS    += -lgfortran

  else
    ifdef BOPT
      PETSC_DIR ?=$(HOME)/lib/petsc-3.0.0-p7

      PETSC_ARCH = linux64_openmpi

      #Determined via modules
      HDF5 = t
      HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(PHDF5_LIB)
      HDF5_INC  = $(PHDF5_INC)

      #Determined via modules
      ADIOS = t
      ADIOS_LIBS = -I$(ADIOS_INC) -L$(ADIOS_LIB) -ladiosf -ladiosreadf -lmxml

    else
      LIBS += -L~/lib -llapack -lblas

      ifeq ($(FC),gfortran)

        #Determined via modules
        HDF5 = t
        HDF5_LIBS =-L$(SHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  = $(SHDF5_LIB)
        HDF5_INC  = $(SHDF5_INC)

        ADIOS = f
      endif
    endif

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

endif

##################
# quevedo (ORNL) #
##################

ifeq ($(HOST),quevedo)

  FC = gfortran

  ifdef SAMR

    CC = mpicc
    CXX = mpicxx
    F77 = gfortran
    override PETSC_DIR  =
#    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/debug
    SAMRSOLVERS_HOME =/opt/samrsolvers-v3/build/optimized
    HDF5 = f
    HDF5_LIBS = 
    LDFLAGS    += -lgfortran

  else
    ifdef BOPT
      PETSC_DIR ?=$(HOME)/lib/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
        PETSC_ARCH = $(PETSC_ARCH_O)
      else
        PETSC_ARCH = $(PETSC_ARCH_g)
      endif

      HDF5 = t

      HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(PHDF5_LIB)
      HDF5_INC  = $(PHDF5_INC)

      ADIOS = t

      #Determined via modules
#      ADIOS_HOME = $(HOME)/lib/adios-1.3
#      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
#      ADIOS_VER  = 1.3
      ADIOS_LIBS = -I$(ADIOS_INC) -L$(ADIOS_LIB) -ladiosf -ladiosreadf -lmxml

    else
      LIBS += -L/usr/lib64 -llapack -lblas
      HDF5 = f
      ifeq ($(FC),f90)
        HDF5 = t

        HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
        HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  = $(HDF5_HOME)/lib
        HDF5_INC  = $(HDF5_HOME)/include

      endif
      ifeq ($(FC),gfortran)
        HDF5 = t

        HDF5_HOME =/usr/lib64
        HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  =$(HDF5_HOME)/gfortran/modules
        HDF5_INC  = /usr/include

        ADIOS = f
      endif
    endif

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

endif

###########
# gongora #
###########

ifeq ($(HOST),gongora)

  FC=gfortran

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = macosx_openmpi_O
    else
      PETSC_ARCH = macosx_openmpi_g
    endif

    HDF5 = t

    HDF5_HOME =/opt/local
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

#    ADIOS = f
    ADIOS = t

    ADIOS_HOME = $(HOME)/lib/adios-1.3.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L$(HOME)/lib/mxml-2.7/lib -lmxml
    ADIOS_VER  = 1.31

  else
    LIBS += -L/usr/lib -llapack -lblas

    HDF5 = t

    HDF5_HOME = $(HOME)/lib/hdf5-1.8.5/gfortran/serial
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/include
    HDF5_INC  =$(HDF5_HOME)/include

    ADIOS = f

  endif

endif

##########
# siesta #
##########

ifeq ($(HOST),siesta)

  FC=gfortran

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi_g
    endif

    HDF5 = t
    HDF5_HOME = /usr/lib64/openmpi
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = /usr/include/openmpi-x86_64
    HDF5_INC  = /usr/include/openmpi-x86_64

    ADIOS = f

    ADIOS_HOME = $(HOME)/lib/adios-1.2.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
    ADIOS_VER  = 1.21

  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f
    ifeq ($(FC),f90)
      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  = $(HDF5_HOME)/include
      HDF5 = t
    endif
    ifeq ($(FC),gfortran)
      HDF5_HOME =/usr/lib64
      HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  =$(HDF5_HOME)/gfortran/modules
      HDF5_INC  = /usr/include
#      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/serial
#      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
#      HDF5_MOD  = $(HDF5_HOME)/lib
#      HDF5_INC  = $(HDF5_HOME)/include
      HDF5 = t
    endif
  endif

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

########
# swim #
########

ifeq ($(HOST),swim)

  FC=gfortran

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = $(PETSC_ARCH_O)
    else
      PETSC_ARCH = $(PETSC_ARCH_g)
    endif

    HDF5 = t

    HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(PHDF5_LIB)
    HDF5_INC  = $(PHDF5_INC)

    ADIOS = t
    ADIOS_LIBS = -I$(ADIOS_INC) -L$(ADIOS_LIB) -ladiosf -ladiosreadf -lmxml

  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f

    HDF5_HOME =/usr/lib64
    HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/gfortran/modules
    HDF5_INC  = /usr/include
    HDF5 = t

    ADIOS = f
  endif

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

##############
# ra22 (RFX) #
##############

ifeq ($(HOST),ra22)

  CPPFLAGS += $(PREPROC)RFX

  FC = pgf95

  override HDF5 = f

  ifdef BOPT
    PETSC_DIR =/ricercatori/ft/petsc-3.0.0-p7
#    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi_g
    endif

    ADIOS = t
    ADIOS_HOME = /ricercatori/ft/adios-1.3
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L/ricercatori/ft/mxml-2.6/lib -lmxml
    ADIOS_VER  = 1.3

  else
    LIBS     += -llapack -lblas

    ADIOS = f
  endif

  NETCDF=f
#  NETCDF_LIBS = -lnetcdff -lnetcdf
#  NETCDF_INC  = -I/usr/include

endif

############
# franklin #
############

ifeq ($(HOST),franklin)
  FC = pgf95
  HDF5_HOME = $(HDF5_PAR_DIR)
#  HDF5_LIBS:= $(HDF5)
  override HDF5 = t
  BOPT = t
endif

#########
# smoky #
#########

ifeq ($(findstring smoky,$(HOST)),smoky)
  FC = pgf95
#  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5 = t
  BOPT = t
endif

##########
# jaguar #
##########

ifeq ($(findstring jaguar,$(HOST)),jaguar)
  FC = pgf95

  BOPT = t
  PETSC_DIR ?= $(CRAY_PETSC_DIR)
  PETSC_ARCH ?= cray-xt

  HDF5 = t
  #HDF5_DIR = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
  HDF5_HOME = $(CRAY_HDF5_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib

  ADIOS = t
  ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) 
  ADIOS_VER  = 1.3
endif

########
# sith #
########

ifeq ($(findstring sith,$(HOST)),sith)
  FC = pgf95

  BOPT = O
  PETSC_DIR ?= /sw/sith/petsc/3.0.0/centos5.5_pgi10.9_opt/petsc-3.0.0-p10
  PETSC_ARCH ?= sith-opt

  HDF5 = t
  #HDF5_DIR = /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2
  HDF5_HOME = $(HDF5_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib -I$(HDF5_DIR)/include

  ADIOS = t
  ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) 
  ADIOS_VER  = 1.3
endif

###########
# chester #
###########

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
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB)
  ADIOS_VER  = 1.3
endif

########
# ewok #
########

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

  ADIOS = t
  #ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_HOME = /ccs/proj/e2e/pnorbert/ADIOS/ewok.noxml
  #ADIOS_LIBS = $(ADIOS_LIB) 
  ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
  ADIOS_VER  = 1.3

  OPTIMIZATION = -O2 -Mbounds
  DEBUG = -g -Mbounds -Mchkptr
  BOPT = t
endif

###########################
# chenhp (Guangye Chen's) #
###########################

ifeq ($(HOST),chenshp)
    LIBS += -L/usr/lib64/atlas -llapack -lblas
    HDF5=f
endif

#dlg-hp overrides 

ifeq ($(HOSTNAME),dlg-hp.ornl.gov)
	PETSC_DIR ?=/home/dg6/code/petsc/petsc-3.1-p7
	HDF5_HOME ?= /usr/lib64/openmpi/lib 
	FC =gfortran
    LIBS += -L ~/code/goto_blas/GotoBLAS2 ~/code/lapack/lapack-3.1.1/lapack_LINUX.a -lgoto2
    HDF5=f
endif

ifeq ($(PIC_CUDA),t) 
	NVCC = nvcc	
	gpu_arch = -arch sm_20
	NVCCFLAGS= --ptxas-options=-v -O2 --compiler-bindir /home/gcc/soft/gcc/gcc-445/bin  #-g -G # 
	CUDA_OPT = -L/usr/local/cuda/lib64 /usr/local/cuda/lib64/libcudart.so  -lcuda -I/usr/local/cuda/include/ 
	CUDA_INCLUDE =  # -I/home/gcc/NVIDIA_GPU_Computing_SDK/C/common/inc
endif

####################
# juelich (JUROPA) #
####################

ifeq ($(findstring jf,$(HOST)),jf)
  FC = ifort

  BOPT = O

  ifdef BOPT
    PETSC_ARCH ?= juropa-basic

    HDF5 = f
    #HDF5_DIR = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
    #HDF5_HOME = $(CRAY_HDF5_DIR)
    #HDF5_LIBS = $(HDF5_FLIB)
    #HDF5_MOD  = $(HDF5_HOME)/lib

    ADIOS = t
    ADIOS_HOME = $(HOME)/lib/adios-1.3
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml -llustreapi
    ADIOS_VER  = 1.3
  else
    HDF5 = f
    ADIOS = f
    LIBS += -llapack -lblas
  endif
endif
