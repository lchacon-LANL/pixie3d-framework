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

HOST ?= $(shell hostname -s)

MACHINE ?=HaveNoIdea

###############################
# lap0086227 (Mark Berrill's) #
###############################

ifeq ($(HOST),lap0086227)
    FC=gfortran
#    FC=vtf77
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

###########
# quevedo #
###########

ifeq ($(HOST),quevedo)
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

###########
# gongora #
###########

ifeq ($(HOST),gongora)
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

##########
# siesta #
##########

ifeq ($(HOST),siesta)
  MACHINE=siesta

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
    HDF5_INC  =-I/usr/include/openmpi-x86_64

    ADIOS = f

    ADIOS_HOME = $(HOME)/lib/adios-1.2.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
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
  NETCDF_INC  = -I/usr/include

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

########
# swim #
########

ifeq ($(HOST),swim)
  MACHINE=swim

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
    HDF5_INC  =-I/usr/include/openmpi-x86_64

    ADIOS = t

    ADIOS_HOME = $(HOME)/lib/adios-1.2.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml

  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f

    HDF5_HOME =/usr/lib64
    HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/gfortran/modules
    HDF5_INC  =-I/usr/include
    HDF5 = t

    ADIOS = f
  endif

  VMEC = t

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

##############
# ra22 (RFX) #
##############

ifeq ($(HOST),ra22)
  MACHINE=ra22

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
    ADIOS_HOME = /ricercatori/ft/adios-1.2.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L/ricercatori/ft/mxml-2.6/lib -lmxml

  else
    LIBS     += -llapack -lblas

    ADIOS = f
  endif

  VMEC = t

  NETCDF=f
#  NETCDF_LIBS = -lnetcdff -lnetcdf
#  NETCDF_INC  = -I/usr/include

endif

############
# franklin #
############

ifeq ($(HOST),franklin)
  MACHINE=franklin

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
  MACHINE=smoky

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

  MACHINE = jaguar
endif

########
# sith #
########

ifeq ($(findstring sith,$(HOST)),sith)
  FC = pgf95

  BOPT = t
  PETSC_DIR ?= /sw/sith/petsc/3.0.0/centos5.5_pgi10.9_opt/petsc-3.0.0-p10
  PETSC_ARCH ?= sith-opt

  HDF5 = t
  #HDF5_DIR = /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2
  HDF5_HOME = $(HDF5_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5_MOD  = $(HDF5_HOME)/lib

  ADIOS = t
  #ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/sith.pgi
  #ADIOS_INC=-I/ccs/proj/e2e/pnorbert/ADIOS/sith.pgi/include -I/sw/sith/mxml/2.6/centos5.5_pgi10.9/include -I/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2/include
  #ADIOS_FLIB=-L/ccs/proj/e2e/pnorbert/ADIOS/sith.pgi/lib -ladiosf -L/sw/sith/mxml/2.6/centos5.5_pgi10.9/lib -lmxml -L/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2/lib -lm -lmxml -lpthread -I/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2/include -L/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2/lib -lhdf5_hl -lhdf5 -L/sw/sith/szip/2.1/centos5.5_pgi10.9/lib -lsz -lz -lm -libverbs

  ADIOS_HOME = $(ADIOS_DIR)
  ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB) 

  MACHINE = sith
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

  MACHINE = chester
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

###########################
# chenhp (Guangye Chen's) #
###########################

ifeq ($(HOST),chenshp)
    MACHINE=chenshp

    LIBS += -L/usr/lib64/atlas -llapack -lblas
    HDF5=f
endif

#dlg-hp overrides 

ifeq ($(HOSTNAME),dlg-hp.ornl.gov)
    MACHINE=dlghp
	PETSC_DIR ?=/home/dg6/code/petsc/petsc-3.1-p7
	HDF5_HOME ?= /usr/lib64/openmpi/lib 
	FC =gfortran
    LIBS += -L ~/code/goto_blas/GotoBLAS2 ~/code/lapack/lapack-3.1.1/lapack_LINUX.a -lgoto2
    HDF5=f
endif

ifeq ($(PIC_CUDA),t) 
	NVCC = nvcc	
	gpu_arch = -arch sm_13
	NVCCFLAGS= --ptxas-options=-v -O2 #-g -G #
	CUDA_OPT = -L/usr/local/cuda/lib64 /usr/local/cuda/lib64/libcudart.so.3.2.16  -lcuda -I/usr/local/cuda/include/ 
	CUDA_INCLUDE = -I/home/gcc/NVIDIA_GPU_Computing_SDK/C/common/inc
endif


