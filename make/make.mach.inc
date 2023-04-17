# Machine configuration

HOST ?= $(shell hostname -s)

ARCH ?= $(shell uname -s)

M_ARCH ?= native

MACHINE ?= $(HOST)

CPPFLAGS += $(PREPROC)$(ARCH)

MPIEXEC ?= mpiexec --oversubscribe

############
# DEFAULTS #
############

FC   ?= gfortran
HDF5 ?= f
ADIOS?= f

###########
# adiosVM #
###########

ifeq ($(HOST),adiosVM)

  override MACHINE = adiosVM

  FC = gfortran

  ifdef BOPT
    PETSC_DIR ?=$(HOME)/Software/petsc-3.0.0-p7

    PETSC_ARCH = linux64_openmpi

    HDF5 = f

    #Determined via modules
    HDF5_LIBS =-L/opt/hdf5-parallel/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm 
    HDF5_MOD  = /opt/hdf5-parallel/lib
    HDF5_INC  = /opt/hdf5-parallel/include

#    ADIOS = f
#    ADIOS_BUFFER_MB = 100
#    ADIOS_METHOD = MPI
#    ADIOS_HOME = /home/adios/Software/adios-1.10.0/mpi-no-hdf5
#    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 
#    ADIOS_VER  = 1.10.0

    ADIOS2 = t
    ADIOS_HOME = /opt/adios2
    ADIOS_LIBS = $(shell $(ADIOS_HOME)/bin/adios2-config -f)
    ADIOS_MOD  = $(ADIOS_HOME)/include/adios2/fortran
    ADIOS_VER  = 2.5.0

  else
    HDF5 = f
    LIBS += -L/home/adios/Software/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netlib-lapack-3.8.0-entm4irbcn2b5io72rmkmivkounzeljd/lib -llapack -lblas

    ADIOS = f
  endif
endif

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
        HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
        HDF5_MOD  = $(HDF5_HOME)/lib
        HDF5_INC  =-I$(HDF5_HOME)/include
    endif
    LIBS += -L/packages/acml-4-4-0/gfortran64/lib -lgfortran
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

  override MACHINE = jetta

  FC ?= gfortran

  OMP_NUM_THREADS = 32

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

  ifeq ($(findstring gfortran,$(FC)),gfortran)
    ifdef BOPT
      PETSC_DIR ?=$(HOME)/lib/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
        PETSC_ARCH = $(PETSC_ARCH_O)
      else
        PETSC_ARCH = $(PETSC_ARCH_g)
      endif

      HDF5 = t

      #Determined via modules
      HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm 
      HDF5_MOD  = $(PHDF5_INC)
      HDF5_INC  = $(PHDF5_INC)

      ifeq ($(PIT),t)
        ADIOS = f      #ADIOS not ready for parallel-in-time yet
      else
        ADIOS = t
      endif

      ADIOS_BUFFER_MB = 100
      ADIOS_METHOD = MPI

      ADIOS_HOME = /nh/u/lchacon/lib/adios-1.10.0
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1
      ADIOS_VER  = 1.10.0

      SPRNG = t

      SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib -lsprng -lstdc++ -lm -lmpi_cxx
      SPRNG_INC  = $(SPRNG_PAR_HOME)/include
    else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = t

      #Determined via modules
      HDF5_LIBS =-L$(SHDF5_LIB) -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
      HDF5_MOD  = $(SHDF5_INC)
      HDF5_INC  = $(SHDF5_INC)

      ADIOS = f

      SPRNG = f

      SPRNG_LIBS = $(SPRNG_SER_HOME)/lib -lsprng -lstdc++ -lm 
      SPRNG_INC = $(SPRNG_SER_HOME)/include
    endif

    NETCDF=t
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

  ifeq ($(findstring ifort,$(FC)),ifort)
    ifdef BOPT
      PETSC_DIR ?=$(HOME)/lib/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
        PETSC_ARCH = $(PETSC_ARCH_O)
      else
        PETSC_ARCH = $(PETSC_ARCH_g)
      endif

      HDF5 = f

      #Determined via modules
      HDF5_LIBS =-L$(PHDF5_LIB) -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz 
      HDF5_MOD  = $(PHDF5_INC)
      HDF5_INC  = $(PHDF5_INC)

      ifeq ($(PIT),t)
        ADIOS = f      #ADIOS not ready for parallel-in-time yet
      else
        ADIOS = f
      endif

      ADIOS_BUFFER_MB = 10
      ADIOS_METHOD = MPI

      ADIOS_HOME = /nh/u/lchacon/lib/adios-1.10.0
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1
      ADIOS_VER  = 1.10.0

      SPRNG = f

      SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib -lsprng -lstdc++ -lmpi_cxx
      SPRNG_INC  = $(SPRNG_PAR_HOME)/include
    else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = f

      #Determined via modules
      HDF5_LIBS =-L$(SHDF5_LIB) -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
      HDF5_MOD  = $(SHDF5_INC)
      HDF5_INC  = $(SHDF5_INC)

      ADIOS = f

      SPRNG = f

      SPRNG_LIBS = $(SPRNG_SER_HOME)/lib -lsprng -lstdc++
      SPRNG_INC = $(SPRNG_SER_HOME)/include
    endif
  endif
  endif
endif

################
# mavis (LANL) #
################

ifeq ($(findstring mavis,$(HOST)),mavis)

  FC = gfortran

  CPPFLAGS += $(PREPROC)sed_nosymlinks

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
      PETSC_DIR ?=/apps/contrib/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
        PETSC_ARCH = linux64_openmpi_gcc-7.2.0
      else
        PETSC_ARCH = linux64_openmpi_gcc-7.2.0_g
      endif

      #Determined via modules
      HDF5 = t
      HDF5_LIBS =-L/apps/hdf5-1.10.2/lib -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
      HDF5_MOD  = /apps/hdf5-1.10.2/include 
      HDF5_INC  = /apps/hdf5-1.10.2/include

      #Determined via modules
      ADIOS = t
      ADIOS_LIBS = -I/apps/contrib/adios-1.10.0/lib/include -L/apps/contrib/adios-1.10.0/lib/lib -ladiosf_v1
      ADIOS_VER  = 1.10.0

    else
      LIBS += -L$(HOME)/lib -llapack -lblas

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

###########
# gongora #
###########

ifeq ($(HOST),gongora)

  FC=gfortran

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

  OMP_NUM_THREADS = 8

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = macosx_openmpi_O
    else
      PETSC_ARCH = macosx_openmpi_g
    endif

    HDF5 = t

#    HDF5_HOME =/opt/local
#    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
#    HDF5_MOD  = $(HDF5_HOME)/include
#    HDF5_INC  = $(HDF5_MOD)
    HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gfortran/parallel
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

    HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gfortran/serial
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

  endif

endif

###########
# quevedo #
###########

ifeq ($(findstring quevedo,$(HOST)),quevedo)

  FC=gfortran

  NETCDF=t
  NETCDF_LIBS = -L/opt/local/lib -lnetcdff -lnetcdf
  NETCDF_INC  = -I/opt/local/include

  OMP_NUM_THREADS = 4

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  ifdef BOPT
    PETSC_DIR =$(HOME)/Software/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = macosx_gcc6.5.0_openmpi3.1.5_O
    else
      PETSC_ARCH = macosx_gcc6.5.0_openmpi3.1.5_g
    endif

    HDF5 = t

    HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gcc-6.5.0/openmpi-3.1.5
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f
    ADIOS_BUFFER_MB = 100
    ADIOS_METHOD = MPI
    ADIOS_HOME = $(HOME)/lib/adios-1.10.0/gcc-6.5.0-openmpi-3.1.5
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1
    ADIOS_VER  = 1.10.0

    ADIOS2 = t
    ADIOS2_METHOD = MPI
    ADIOS2_HOME = $(HOME)/lib/adios-2.5.0/gcc-6.5.0_openmpi-3.1.5
    ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
    ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags)
    ADIOS2_VER  = 2.5.0

    SPRNG = t
    SPRNG_HOME = $(HOME)/lib/sprng5/gcc-6.5.0_openmpi-3.1.5
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC  = $(SPRNG_HOME)/include
      
#    LIBS += -L/opt/local/lib -latlas
  else
    LIBS += -L/opt/local/lib/lapack -llapack -lblas

    HDF5 = t

    HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gcc-6.5.0/serial
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

    ADIOS2 = f

    SPRNG = t
    SPRNG_HOME = $(HOME)/lib/sprng5/gcc-6.5.0_serial
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC  = $(SPRNG_HOME)/include

  endif

endif

#############
# pn2033419 #
#############

ifeq ($(findstring pn2033419,$(HOST)),pn2033419)

  FC=gfortran

  NETCDF=t
  NETCDF_LIBS = -L/opt/local/lib -lnetcdff -lnetcdf
  NETCDF_INC  = -I/opt/local/include

  OMP_NUM_THREADS = 8

  override MACHINE = pn2033419

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  ifdef BOPT
    ifdef MPICH
      MPIEXEC = mpiexec
      
      PETSC_DIR =$(HOME)/Software/petsc-3.17.3

#      CPPFLAGS+= $(PREPROC)USE_MPICH

      ifeq ($(BOPT),O)
        PETSC_ARCH = darwin_gcc-9.4_mpich-3.4.2
      else
        PETSC_ARCH = darwin_gcc-9.4_mpich-3.4.2_g
      endif

      HDF5 = f
      HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gcc-9.4.0-openmpi-3.1.5
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/include
      HDF5_INC  = $(HDF5_MOD)

      ADIOS2 = t
      ADIOS2_METHOD = MPI
      ADIOS2_HOME = $(HOME)/lib/adios-2.8.3/gcc-9.4.0_mpich-3.4.2
      ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
      ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags)
      ADIOS2_VER  = 2.8.3

      SPRNG = t
      SPRNG_HOME = $(HOME)/lib/sprng5/gcc-9.4.0_mpich-3.4.2
      SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
      SPRNG_INC  = $(SPRNG_HOME)/include
    else
      PETSC_DIR =$(HOME)/Software/petsc-3.17.3

      ifeq ($(BOPT),O)
        PETSC_ARCH = darwin_gcc-9.4_openmpi-3.1.5
      else
        PETSC_ARCH = darwin_gcc-9.4_openmpi-3.1.5_g
      endif

      HDF5 = t
      HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gcc-9.4.0-openmpi-3.1.5
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/include
      HDF5_INC  = $(HDF5_MOD)

      ADIOS2 = t
      ADIOS2_METHOD = MPI
#      ADIOS2_HOME = $(HOME)/lib/adios-2.5.0/gcc-9.4.0_openmpi-3.1.5
      ADIOS2_HOME = $(HOME)/lib/adios-2.8.3/gcc-9.4.0_openmpi-3.1.5
      ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
      ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags)
#      ADIOS2_VER  = 2.5.0
      ADIOS2_VER  = 2.8.3

      SPRNG = t
      SPRNG_HOME = $(HOME)/lib/sprng5/gcc-9.4.0_openmpi-3.1.5
      SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
      SPRNG_INC  = $(SPRNG_HOME)/include
    endif
      
  else
    LIBS += -llapack -lblas

    HDF5 = t

    HDF5_HOME = $(HOME)/lib/hdf5-1.8.13/gcc-9.4.0
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

    ADIOS2 = f

    SPRNG = t
    SPRNG_HOME = $(HOME)/lib/sprng5/gcc-9.4.0
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC  = $(SPRNG_HOME)/include

  endif

endif
##############
# quadrivium #
##############

ifeq ($(HOST),quadrivium)

  FC=gfortran

  PIT=t
  VECPOT=t

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi_g
    endif

    HDF5 = t

    HDF5_HOME = $(HOME)/local
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

    ADIOS_HOME = $(HOME)/lib/adios-1.3.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L$(HOME)/lib/mxml-2.7/lib -lmxml
    ADIOS_VER  = 1.31

  else

    override LDFLAGS=

    LIBS += -L/state/partition1/apps/lapack/lib -llapack -lrefblas

    HDF5 = t

    HDF5_HOME = $(HOME)/local2
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/include
    HDF5_INC  =$(HDF5_HOME)/include

    ADIOS = f

  endif

endif

#######################
# Maria's Macbook Air #
#######################

ifeq ($(findstring Marias,$(HOST)),Marias)

  FC=gfortran

  MACHINE=Macbook-Air

  NETCDF=f
  NETCDF_LIBS = -L/opt/local/lib -lnetcdff -lnetcdf
  NETCDF_INC  = -I/opt/local/include

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  WORKDIR = $(HOME)/work

  OMP_NUM_THREADS = 4

  ifdef BOPT
    PETSC_DIR =$(WORKDIR)/lib/petsc-3.17.3

    ifeq ($(BOPT),O)
      PETSC_ARCH = darwin_gcc-12.2.0_openmpi-3.1.6
    else
      PETSC_ARCH = darwin_gcc-12.2.0_openmpi-3.1.6_g
    endif

    HDF5 = t

    HDF5_HOME = $(WORKDIR)/lib/hdf5-1.8.13/gcc-12.2.0_openmpi-3.1.6
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

#    ADIOS = t
#    ADIOS_HOME = $(WORKDIR)/lib/adios-2.8.2/gcc-12.1.0_openmpi-3.1.6
#    ADIOS_LIBS = -L$(ADIOS_HOME)/lib -ladiosf_v1
#    ADIOS_INC  = $(ADIOS_HOME)/include
#    ADIOS_VER  = 1.10.0

    ADIOS2 = t
    ADIOS2_METHOD = MPI
    ADIOS2_HOME = $(WORKDIR)/lib/adios-2.5.0/gcc-12.2.0_openmpi-3.1.6
    ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
    ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
    ADIOS2_VER  = 2.5.0
     
    LIBS += -L/opt/local/lib -latlas
  else
    LIBS += -L/opt/local/lib -llapack -lblas -latlas

    HDF5 = t

    HDF5_HOME = $(WORKDIR)/lib/hdf5-1.8.13/gcc12.1
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/include
    HDF5_INC  =$(HDF5_HOME)/include

    ADIOS = f

  endif

endif

##############
# ra22 (RFX) #
##############

ifeq ($(HOST),ra22)

  CPPFLAGS += $(PREPROC)RFX $(PREPROC)sed_nosymlinks

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
   FC = gfortran

# BOPT = O
#PETSC_DIR = /sw/sith/petsc/3.0.0/rhel6_pgi13.4_debug/petsc-3.0.0-p10
 #/sw/sith/petsc/3.0.0/centos5.5_pgi10.9_opt/petsc-3.0.0-p10
  PETSC_DIR = /sw/sith/petsc/3.0.0/rhel6_gnu4.8.2_opt/petsc-3.0.0-p10
   PETSC_ARCH = sith-opt

   HDF5 = t
 #HDF5_DIR = /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.2
#HDF5_DIR = /sw/sith/hdf5/1.8.10/rhel6_pgi12.8
#HDF5_DIR = /sw/sith/hdf5/1.8.10/rhel6_gnu4.8.2_ompi1.8.3
   HDF5_HOME = $(HDF5_DIR)
   HDF5_LIBS = $(HDF5_FLIB)
   HDF5_MOD  = $(HDF5_HOME)/lib #-I$(HDF5_DIR)/include
   HDF5_INC  = $(HDF5_DIR)/include

   ADIOS = t
#ADIOS_HOME = $(ADIOS_DIR)
#  ADIOS_INC=-I/lustre/atlas2/csc143/proj-shared/tian/adios-1.3.1-install/include -I/sw/redhat6/mxml/2.9/rhel6_gnu4.7.1/includ
#   ADIOS_FLIB=-L/lustre/atlas2/csc143/proj-shared/tian/adios-1.3.1-install/lib -ladiosf -L/sw/redhat6/mxml/2.9/rhel6_gnu4.7.1/lib -lmxml -lm -lmxml -lpthread -libverbs

  ADIOS_HOME = /lustre/atlas2/csc143/proj-shared/tian/adios-1.10.0-install
  ADIOS_INC=-I/lustre/atlas2/csc143/proj-shared/tian/adios-1.10.0-install/include -I/sw/redhat6/dataspaces/1.6.0/rhel6.6_gnu4.8.2/include -DZLIB -DBZIP2 -I/sw/sith/bzip2/1.0.6/rhel6_gnu4.7.1/include
  ADIOS_FLIB=-L/lustre/atlas2/csc143/proj-shared/tian/adios-1.10.0-install/lib -ladiosf_v1 -L/sw/redhat6/dataspaces/1.6.0/rhel6.6_gnu4.8.2/lib -L/sw/sith/bzip2/1.0.6/rhel6_gnu4.7.1/lib64 -ldspaces -ldscommon -ldart -lrdmacm -libverbs -llustreapi -lz -lbz2 -lm -lrt
   ADIOS_LIBS = $(ADIOS_INC) $(ADIOS_FLIB)
   ADIOS_VER  = 1.10
endif

##########################
# Guangye Chen's desktop #
##########################

ifeq ($(HOST),t1243055)

  FC=gfortran

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  WORKDIR = $(HOME)

  ifdef BOPT
    PETSC_DIR =$(WORKDIR)/soft/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi
    endif

    HDF5 = f

    HDF5_HOME =/opt/local
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

    ADIOS_HOME = $(WORKDIR)/lib/adios-1.3.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L$(WORKDIR)/lib/mxml-2.7/lib -lmxml
    ADIOS_VER  = 1.31

    SPRNG = t

    SPRNG_HOME = /home/gchen/Research/example/sprng/sprng5/parallel
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC  = $(SPRNG_HOME)/include
  else
    LIBS += -llapack -lblas

    HDF5 = f

    HDF5_HOME = $(WORKDIR)/lib/hdf5-1.8.7/gfortran/serial
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/include
    HDF5_INC  =$(HDF5_HOME)/include

    ADIOS = f

    SPRNG = t

    SPRNG_HOME = /home/gchen/Research/example/sprng/sprng5/serial
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC = $(SPRNG_HOME)/include
  endif

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
  OPTIMIZATION = -O2

#  BOPT = O

  CPPFLAGS += $(PREPROC)RFX

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

####################
# helios #
####################

ifeq ($(findstring helios,$(HOST)),helios)
  FC = ifort
  OPTIMIZATION = -O3 -xavx

  CPPFLAGS += $(PREPROC)RFX

  ifdef BOPT
#    PETSC_DIR =/csc/project/ENEATEST/PIXIE3D/petsc-3.0.0-p7
    PETSC_DIR =/csc/project/PIXIE3D/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_intelmpi
    else
      PETSC_ARCH = linux64_intelmpi_g
    endif

    HDF5 = f
    #HDF5_DIR = /sw/xt5/hdf5/1.6.7/cnl2.1_pgi7.2.3_par
    #HDF5_HOME = $(CRAY_HDF5_DIR)
    #HDF5_LIBS = $(HDF5_FLIB)
    #HDF5_MOD  = $(HDF5_HOME)/lib

    ADIOS = t
    ADIOS_HOME = /csc/softs/nccs/adios-1.3.1/intel-12.0.5.220/intelmpi-4.0.3/default/
    #ADIOS_HOME = $(ADIOS_DIR)
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml -llustreapi
    ADIOS_VER  = 1.31
    #ADIOS_VER  = 1.50 # error at link time
    #ADIOS_VER  = 1.70 # error at link time
  else
    HDF5 = f
    ADIOS = f
#    LIBS += -llapack -lblas
    LIBS += -mkl
  endif

  override MACHINE = helios
endif

################
# Will Taitano's Desktop(LANL) #
################

ifeq ($(HOST),t1243000)

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
        PETSC_ARCH = linux64_openmpi
      else
        PETSC_ARCH = linux64_openmpi_g
      endif

      HDF5 = t

      #Determined via modules
      HDF5_LIBS =-L/usr/lib64/openmpi/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = /usr/include/openmpi-x86_64
      HDF5_INC  = /usr/include/openmpi-x86_64

      ADIOS = t

      ADIOS_HOME = $(HOME)/lib/adios-1.31
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -lmxml
      ADIOS_VER  = 1.31

      #Determined via modules
      HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = /usr/local/gfortran/
      HDF5_INC  = /usr/include

    else
      HDF5 = f
      ADIOS = f
      LIBS += -llapack -lblas
    endif

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

#    ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
  endif

endif

##################
# Wolf turquoise #
##################

ifeq ($(findstring wf,$(HOST)),wf)

  override MACHINE = wolf

  ifeq ($(findstring gfortran,$(FC)),gfortran)

    NETCDF=t
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

    ifdef BOPT
     PETSC_DIR =/turquoise/usr/projects/magreconn/soft/petsc-3.0.0-p7

     ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_wolf_openmpi-1.10.5_gcc-5.3.0
     else
      PETSC_ARCH = linux64_wolf_openmpi-1.10.5_gcc-5.3.0_dgb
     endif

     HDF5 = t

     HDF5_HOME = /usr/projects/hpcsoft/toss2/wolf/hdf5/1.8.16_gcc-5.3.0_openmpi-1.10.5
     HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_HOME)/include
     HDF5_INC  = $(HDF5_MOD)

     ADIOS = t

     ADIOS_BUFFER_MB = 1000
     ADIOS_METHOD = MPI_LUSTRE

     ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-wolf-gcc-5.3.0_openmpi-1.10.5
     ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
     ADIOS_VER  = 1.10.0

     SPRNG = t

     SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/wolf-gcc5.3.0-openmpi1.10.5/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_PAR_HOME)/include
    else
     HDF5 = f
     ADIOS = f
     LIBS += -llapack -lblas
    endif
  endif

  ifeq ($(findstring ifort,$(FC)),ifort)

    NETCDF=f

    ifdef BOPT
       PETSC_DIR =/turquoise/usr/projects/magreconn/soft/petsc-3.0.0-p7

       ifeq ($(BOPT),O)
        PETSC_ARCH = linux64_wolf_intel17.0.1_openmpi1.10.5
       else
        PETSC_ARCH = linux64_wolf_intel17.0.1_openmpi1.10.5_dbg
       endif

       HDF5 = t

       HDF5_HOME = /usr/projects/hpcsoft/toss2/wolf/hdf5/1.8.16_intel-17.0.1_openmpi-1.10.5
       HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
       HDF5_MOD  = $(HDF5_HOME)/include
       HDF5_INC  = $(HDF5_MOD)

       ADIOS = t

       ADIOS_BUFFER_MB = 1000
       ADIOS_METHOD = MPI_LUSTRE

       ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-wolf-intel-17.0.1_openmpi-1.10.5
       ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
       ADIOS_VER  = 1.10.0

       SPRNG = t

       SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
       SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/wolf-intel17.0.1-openmpi1.10.5/lib -lsprng -lstdc++ -lm
       SPRNG_INC  = $(SPRNG_PAR_HOME)/include

    else
      HDF5 = f
      ADIOS = f
      LIBS += -llapack -lblas
    endif

  endif

endif

##########
# Darwin #
##########

ifeq ($(findstring darwin,$(HOST)),darwin)

  override MACHINE = darwin

  ifeq ($(findstring gfortran,$(FC)),gfortran)

    NETCDF=f
    NETCDF_LIBS = -lnetcdff -lnetcdf
    NETCDF_INC  = -I/usr/include

    ifdef BOPT
     PETSC_DIR =$(HOME)/software/petsc-3.0.0-p7

     ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi3.1.0_gcc7.3.0
     else
      PETSC_ARCH = linux64_openmpi3.1.0_gcc7.3.0_g
     endif

     HDF5 = f

     HDF5_HOME = /usr/projects/hpcsoft/toss2/wolf/hdf5/1.8.16_gcc-5.3.0_openmpi-1.10.5
     HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_HOME)/include
     HDF5_INC  = $(HDF5_MOD)

     ADIOS = f

     ADIOS_BUFFER_MB = 1000
     ADIOS_METHOD = MPI_LUSTRE

     ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-wolf-gcc-5.3.0_openmpi-1.10.5
     ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
     ADIOS_VER  = 1.10.0

     SPRNG = f

     SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/wolf-gcc5.3.0-openmpi1.10.5/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_PAR_HOME)/include
    else
     HDF5 = f
     ADIOS = f
     LIBS += -llapack -lblas
    endif
  endif

  ifeq ($(findstring ifort,$(FC)),ifort)

    NETCDF=f

    ifdef BOPT
       PETSC_DIR =$(HOME)/software/petsc-3.0.0-p7

       ifeq ($(BOPT),O)
        PETSC_ARCH = linux64_openmpi3.1.0_intel17.0.6
       else
        PETSC_ARCH = linux64_openmpi3.1.0_intel17.0.6_g
       endif

       HDF5 = f

       HDF5_HOME = /usr/projects/hpcsoft/toss2/wolf/hdf5/1.8.16_intel-17.0.1_openmpi-1.10.5
       HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
       HDF5_MOD  = $(HDF5_HOME)/include
       HDF5_INC  = $(HDF5_MOD)

       ADIOS = f

       ADIOS_BUFFER_MB = 1000
       ADIOS_METHOD = MPI

       ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-wolf-intel-17.0.1_openmpi-1.10.5
       ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
       ADIOS_VER  = 1.10.0

       SPRNG = f

       SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
       SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/wolf-intel17.0.1-openmpi1.10.5/lib -lsprng -lstdc++ -lm
       SPRNG_INC  = $(SPRNG_PAR_HOME)/include

    else
      HDF5 = f
      ADIOS = f
      LIBS += -llapack -lblas
    endif

  endif

endif

#####################
# Badger turquoise #
#####################

ifeq ($(findstring ba,$(HOST)),ba)

  override MACHINE = badger

  ifeq ($(findstring gfortran,$(FC)),gfortran)

    NETCDF=f

    ifdef BOPT
     PETSC_DIR =/turquoise/usr/projects/magreconn/soft/petsc-3.17.3

     ifdef MPICH_FC
       MPIEXEC = srun

       CPPFLAGS+= $(PREPROC)USE_MPICH
       ifeq ($(BOPT),O)
         PETSC_ARCH = badger_gcc-9.3.0_mvapich2-2.3.6
       else
         PETSC_ARCH = badger_gcc-9.3.0_mvapich2-2.3.6_g
       endif

       ADIOS = f
       ADIOS2 = f
       SPRNG = f
     else
       ifeq ($(BOPT),O)
         PETSC_ARCH = badger_gcc-9.3.0_openmpi-3.1.6
       else
         PETSC_ARCH = badger_gcc-9.3.0_openmpi-3.1.6
       endif

       ADIOS = f
       ADIOS_BUFFER_MB = 1000
       ADIOS_METHOD = MPI_LUSTRE
       ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-badger-gcc-9.3.0_openmpi-3.1.6
       ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
       ADIOS_VER  = 1.10.0

       ADIOS2 = t
       ADIOS2_METHOD = MPI_LUSTRE
       ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios-2.5/lib/badger_gcc-9.3.0_openmpi-3.1.6
#       ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios2-latest/lib/badger_gcc-9.3.0_openmpi-3.1.6
       ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
       ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
       ADIOS2_VER  = 2.5.0
#       ADIOS2_VER  = 2.7.1
     
       SPRNG = t

       SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
       SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/badger_gcc-9.3.0_openmpi-3.1.6/lib -lsprng -lstdc++ -lm
       SPRNG_INC  = $(SPRNG_PAR_HOME)/include
     endif

     HDF5 = t

     HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_ROOT)/include
     HDF5_INC  = $(HDF5_MOD)

    else
     LIBS += -L/usr/lib -llapack -lblas

     HDF5 = f

     ADIOS = f

     SPRNG = t

     SPRNG_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_HOME)/lib/badger_gcc-9.3.0/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_HOME)/include
    endif

  endif

  ifeq ($(findstring ifort,$(FC)),ifort)

    NETCDF=f

    ifeq ($(MPI_NAME),openmpi)
     ifdef BOPT
      PETSC_DIR =/turquoise/usr/projects/magreconn/soft/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
       PETSC_ARCH = linux64_badger_openmpi-3.1.6_intel-19.0.4
      else
       PETSC_ARCH = linux64_badger_openmpi-3.1.6_intel-19.0.4_g
      endif

      HDF5 = t

      HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_ROOT)/include
      HDF5_INC  = $(HDF5_MOD)

      ADIOS = f
      ADIOS_BUFFER_MB = 1000
      ADIOS_METHOD = MPI_LUSTRE
      ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-badger-intel-19.0.4_openmpi-3.1.6
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
      ADIOS_VER  = 1.10.0

      ADIOS2 = t
      ADIOS2_METHOD = MPI_LUSTRE
      ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios-2.5/lib/badger_intel-19.0.4_openmpi-3.1.6
#      ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios2-latest/lib/badger_intel-19.0.4_openmpi-3.1.6
      ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
      ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
      ADIOS2_VER  = 2.5.0
#      ADIOS2_VER  = 2.7.1
     
      SPRNG = t

      SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
      SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/badger_intel-19.0.4_openmpi-3.1.6/lib -lsprng -lstdc++
      SPRNG_INC  = $(SPRNG_PAR_HOME)/include
     else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = f

      ADIOS = f

      SPRNG = t

      SPRNG_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
      SPRNG_LIBS = $(SPRNG_HOME)/lib/badger_intel-19.0.4/lib -lsprng -lstdc++
      SPRNG_INC  = $(SPRNG_HOME)/include

     endif

    endif

    ifeq ($(MPI_NAME),intel-mpi)
     ifdef BOPT
      PETSC_DIR =/turquoise/usr/projects/magreconn/soft/petsc-3.0.0-p7

      ifeq ($(BOPT),O)
       PETSC_ARCH = linux64_grizzly_intel17.0.1_intelmpi2017.1
      else
       PETSC_ARCH = linux64_grizzly_intel17.0.1_intelmpi2017.1_dbg
      endif

      HDF5 = f

#      HDF5_ROOT = /usr/projects/hpcsoft/toss3/grizzly/hdf5/1.8.16_gcc-5.3.0_openmpi-1.10.3
#      HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
#      HDF5_MOD  = $(HDF5_ROOT)/include
#      HDF5_INC  = $(HDF5_MOD)

      ADIOS = f

      ADIOS_BUFFER_MB = 1000
      ADIOS_METHOD = MPI_LUSTRE

      ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-grizzly-intel-17.0.1_intelmpi-2017.1
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
      ADIOS_VER  = 1.10.0
     endif
   endif

  endif

endif

#####################
# Chicoma turquoise #
#####################

ifeq ($(findstring ch,$(HOST)),ch)

  override MACHINE = chicoma
  override FC = ftn
  
  MPIEXEC = srun

  CPPFLAGS+= $(PREPROC)USE_MPICH

  ifeq ($(findstring gfortran,$(shell $(FC) 2>&1)),gfortran)

    NETCDF=f

    ifdef BOPT
     PETSC_DIR =/usr/projects/magreconn/soft/petsc-3.17.3

     ifeq ($(BOPT),O)
      PETSC_ARCH = chicoma_gcc-12.1.0_mpich-8.1.21
     else
      PETSC_ARCH = chicoma_gcc-12.1.0_mpich-8.1.21_g
     endif

     HDF5 = t

     HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_ROOT)/include
     HDF5_INC  = $(HDF5_MOD)

     ADIOS2 = t
     ADIOS2_METHOD = MPI_LUSTRE
     ADIOS2_HOME = /usr/projects/magreconn/soft/adios-2.5/lib/chicoma_gnu-12.1.0_mpich-8.1.21
     ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
     ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
     ADIOS2_VER  = 2.5.0

     SPRNG = t

     SPRNG_PAR_HOME = /usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/chicoma_gcc-12.1.0_mpich-8.1.21/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_PAR_HOME)/lib/chicoma_gcc-12.1.0_mpich-8.1.21/include

    else
      HDF5 = f
      ADIOS = f
      LIBS += -llapack -lblas
    endif

  endif

  ifeq ($(findstring ifort,$(shell $(FC) 2>&1)),ifort)

    NETCDF=f

    ifdef BOPT
     PETSC_DIR =/usr/projects/magreconn/soft/petsc-3.17.3

     ifeq ($(BOPT),O)
      PETSC_ARCH = chicoma_intel-21.5.0_mpich-8.1.14
     else
      PETSC_ARCH = chicoma_intel-21.5.0_mpich-8.1.14_g
     endif

     HDF5 = t

     HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_ROOT)/include
     HDF5_INC  = $(HDF5_MOD)

     ADIOS2 = f
     ADIOS2_METHOD = MPI_LUSTRE
     ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios-2.5/lib/badger_intel-19.0.4_openmpi-3.1.6
     ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
     ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
     ADIOS2_VER  = 2.5.0

     SPRNG = t

     SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/chicoma_intel-21.0.2_mpich-8.1.14/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_PAR_HOME)/include

    else
      HDF5 = f
      ADIOS = f
      LIBS += -llapack -lblas
    endif

  endif
  
  ifeq ($(findstring ftn,$(shell $(FC) 2>&1)),ftn)

    NETCDF=f

    ifdef BOPT
     PETSC_DIR =/usr/projects/magreconn/soft/petsc-3.17.3

     ifeq ($(BOPT),O)
      PETSC_ARCH = chicoma_cray-13.0.2_mpich-8.1.14
     else
      PETSC_ARCH = chicoma_cray-13.0.2_mpich-8.1.14_g
     endif

     HDF5 = t

     HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm
     HDF5_MOD  = $(HDF5_ROOT)/include
     HDF5_INC  = $(HDF5_MOD)

     ADIOS2 = f
     ADIOS2_METHOD = MPI_LUSTRE
     ADIOS2_HOME = /turquoise/usr/projects/magreconn/soft/adios-2.5/lib/badger_intel-19.0.4_openmpi-3.1.6
     ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
     ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags) 
     ADIOS2_VER  = 2.5.0

     SPRNG = f

     SPRNG_PAR_HOME = /turquoise/usr/projects/magreconn/soft/sprng5
     SPRNG_LIBS = $(SPRNG_PAR_HOME)/lib/grizzly-gcc5.3.0-openmpi1.10.5/lib -lsprng -lstdc++ -lm
     SPRNG_INC  = $(SPRNG_PAR_HOME)/include

    endif

  endif
endif

########
# Snow #
########

ifeq ($(findstring sn,$(HOST)),sn)

  override MACHINE = snow

  ifeq ($(findstring ifort,$(shell $(FC) 2>&1)),ifort)

    NETCDF=f

    ifeq ($(MPI_NAME),openmpi)
     ifdef BOPT
      PETSC_DIR =/usr/projects/DPIC/Software/petsc-3.17.3

      ifeq ($(BOPT),O)
       PETSC_ARCH = snow_intel-21.2.0_openmpi-4.1.1
      else
       PETSC_ARCH = snow_intel-21.2.0_openmpi-4.1.1_g
      endif

      HDF5 = f

      HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz
      HDF5_MOD  = $(HDF5_ROOT)/include
      HDF5_INC  = $(HDF5_MOD)

      ADIOS = f
      ADIOS_BUFFER_MB = 1000
      ADIOS_METHOD = MPI_LUSTRE
      ADIOS_HOME = /turquoise/usr/projects/magreconn/soft/adios-1.10.0/lib-grizzly-intel-17.0.1_openmpi-1.10.5
      ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -llustreapi
      ADIOS_VER  = 1.10.0

      ADIOS2 = t
      ADIOS2_METHOD = MPI_LUSTRE
      ADIOS2_HOME = /usr/projects/DPIC/Software/adios-2.5/lib/snow_intel-21.2.0_openmpi-4.1.1
      ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
      ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags)
      ADIOS2_VER  = 2.5.0
      
      SPRNG = t

      SPRNG_PAR_HOME = /usr/projects/DPIC/Software/sprng5
      SPRNG_LIBS = $(SPRNG_PAR_HOME)/snow_intel-21.2.0-openmpi-4.1.1/lib -lsprng -lstdc++
      SPRNG_INC  = $(SPRNG_PAR_HOME)/include
     else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = f

      ADIOS = f

      SPRNG = f

     endif

    endif

  endif

  ifeq ($(findstring gfortran,$(shell $(FC) 2>&1)),gfortran)

    NETCDF=f

    ifeq ($(MPI_NAME),openmpi)
     ifdef BOPT
      PETSC_DIR =/usr/projects/DPIC/Software/petsc-3.17.3

      ifeq ($(BOPT),O)
       PETSC_ARCH = snow_gcc-10.3.0_openmpi-4.1.1
      else
       PETSC_ARCH = snow_gcc-10.3.0_openmpi-4.1.1_g
      endif

      HDF5 = f

      HDF5_LIBS =-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz
      HDF5_MOD  = $(HDF5_ROOT)/include
      HDF5_INC  = $(HDF5_MOD)

      ADIOS = f

      ADIOS2 = t
      ADIOS2_METHOD = MPI_LUSTRE
      ADIOS2_HOME = /usr/projects/DPIC/Software/adios-2.5/lib/gcc-10.3.0_openmpi-4.1.1
      ADIOS2_LIBS = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-libs)
      ADIOS2_MOD  = $(shell $(ADIOS2_HOME)/bin/adios2-config --fortran-flags)
      ADIOS2_VER  = 2.5.0

      SPRNG = t

      SPRNG_PAR_HOME = /usr/projects/DPIC/Software/sprng5
      SPRNG_LIBS = $(SPRNG_PAR_HOME)/snow_gcc-10.3.0_openmpi-4.1.1/lib -lsprng -lstdc++
      SPRNG_INC  = $(SPRNG_PAR_HOME)/include
     else
      LIBS += -L/usr/lib -llapack -lblas

      HDF5 = f

      ADIOS = f

      SPRNG = f

     endif

    endif

  endif

endif

####################
# marconi #
####################

ifeq ($(findstring marconi_work,$(WORK)),marconi_work)

#  FC = ifort
 FC = gfortran
 OPTIMIZATION = -O3

#  CPPFLAGS += $(PREPROC)RFX
# togliere per casi VMEC così pixplot produce sia cnv sia cov

 ifdef BOPT
   PETSC_DIR = /marconi_work/FUA33_PIXIE3D/petsc-3.0.0-p7

   ifeq ($(BOPT),O)
# ifort/intelmpi
#      PETSC_ARCH = linux64_intelmpi
# gfortran/openmpi
     PETSC_ARCH = linux64_openmpi_gnu
# ifort/openmpi
#      PETSC_ARCH = linux64_openmpi_ifort
   else
     PETSC_ARCH = linux64_intelmpi_g
   endif

   HDF5 = f

#   HDF5 = t
#   HDF5_HOME = /cineca/prod/opt/libraries/hdf5/1.8.17/intelmpi--2017--binary
#   HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
#   HDF5_MOD  = $(HDF5_HOME)/include
#   HDF5_INC  = $(HDF5_MOD)

#   HDF5 = t
#   HDF5_HOME = /cineca/prod/opt/libraries/hdf5/1.8.18/openmpi--1-10.3--gnu--6.1.0
#   HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
#   HDF5_MOD  = $(HDF5_HOME)/include
#   HDF5_INC  = $(HDF5_MOD)

   ADIOS = f

   #ADIOS = t
   #ADIOS_HOME = /marconi_work/FUA11_PIXIE3D/adios-1.3.1/linux64_intelmpi
   #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L/marconi_work/FUA11_PIXIE3D/mxml-2.10/linux64_intelmpi/lib -lmxml
   #ADIOS_VER  = 1.31

   #ADIOS = t
   #ADIOS_HOME = /marconi_work/FUA11_PIXIE3D/adios-1.10.0/linux64_intelmpi
   #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -L/marconi_work/FUA11_PIXIE3D/mxml-2.10/linux64_intelmpi/lib -lmxml
   #ADIOS_VER  = 1.10

   #ADIOS = t
   #ADIOS_HOME = /marconi_work/FUA11_PIXIE3D/adios-1.3.1/linux64_openmpi_gnu
   #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L/marconi_work/FUA11_PIXIE3D/mxml-2.10/linux64_openmpi_gnu/lib -lmxml
   #ADIOS_VER  = 1.31


   #ADIOS = t
   #ADIOS_HOME = /marconi_work/FUA11_PIXIE3D/adios-1.10.0/linux64_openmpi_gnu
   #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -L/marconi_work/FUA11_PIXIE3D/mxml-2.10/linux64_openmpi_gnu/lib -lmxml
   #ADIOS_VER  = 1.10

   #ADIOS = t
   #ADIOS_HOME = /marconi_work/FUA11_PIXIE3D/adios-1.10.0/linux64_openmpi_ifort
   #ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf_v1 -L/marconi_work/FUA11_PIXIE3D/mxml-2.10/linux64_openmpi_ifort/lib -lmxml
   #ADIOS_VER  = 1.10

   NETCDF=t
   NETCDF_LIBS = -lnetcdff -lnetcdf
   NETCDF_INC  = -I/usr/include

 else
   HDF5 = f
   ADIOS = f
   LIBS += -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas
 endif

   NETCDF=f

#   NETCDF=t
#   NETCDF_LIBS = -L/cineca/prod/opt/libraries/netcdf/4.4.1/gnu--6.1.0/lib -lnetcdf -L/cineca/prod/opt/libraries/netcdff/4.4.4/gnu--6.1.0/lib -lnetcdff
#   NETCDF_INC = -I/cineca/prod/opt/libraries/netcdf/4.4.1/gnu--6.1.0/include -I/cineca/prod/opt/libraries/netcdff/4.4.4/gnu--6.1.0/include

#   NETCDF=t
#    NETCDF_LIBS = -L/cineca/prod/opt/libraries/netcdf/4.4.1/intel--pe-xe-2017--binary/lib -lnetcdff -lnetcdf
#    NETCDF_INC = -I/cineca/prod/opt/libraries/netcdf/4.4.1/intel--pe-xe-2017--binary/include

 MACHINE = marconi
endif

###############################
#  (Koshkarov O.  ) #
###############################

ifeq ($(HOST),cascabel)
    FC=gfortran
    #FC=mpif90
    #CC=mpicc
    #CXX=mpic++
    HDF5 = f
    HDF5_HOME = /home/koshmar/.local/hdf5
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  = $(HDF5_HOME)/include
    LIBS +=  -llapack -lblas
endif

ifeq ($(HOST),kigal)
    FC=gfortran
    #FC=mpif90
    #CC=mpicc
    #CXX=mpic++
    HDF5 = t
    HDF5_HOME = /home/koshmar/.local/hdf5
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  = $(HDF5_HOME)/include
    LIBS +=  -llapack -lblas
    
    # What is CPPFLAGS? is it FFLAGS?
endif

###############################
#  (drava LANL) #
###############################
ifeq ($(HOST),drava)

  FC=gfortran

  NETCDF=f
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)

  WORKDIR = $(HOME)

  ifdef BOPT
    PETSC_DIR =$(WORKDIR)/soft/petsc-3.0.0-p7

    ifeq ($(BOPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi
    endif

    HDF5 = f

    HDF5_HOME =/opt/local
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/include
    HDF5_INC  = $(HDF5_MOD)

    ADIOS = f

    ADIOS_HOME = $(WORKDIR)/lib/adios-1.3.1
    ADIOS_LIBS = -I$(ADIOS_HOME)/include -L$(ADIOS_HOME)/lib -ladiosf -ladiosreadf -L$(WORKDIR)/lib/mxml-2.7/lib -lmxml
    ADIOS_VER  = 1.31

    SPRNG = t

    SPRNG_HOME = /home/scratch/guangyechen/t1243055/example/sprng/sprng5/parallel
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC  = $(SPRNG_HOME)/include
  else
    LIBS += -llapack -lblas

    HDF5 = f

    HDF5_HOME = $(WORKDIR)/lib/hdf5-1.8.7/gfortran/serial
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  =$(HDF5_HOME)/include
    HDF5_INC  =$(HDF5_HOME)/include

    ADIOS = f

    SPRNG = t

    SPRNG_HOME = /home/scratch/guangyechen/t1243055/example/sprng/sprng5/serial
    SPRNG_LIBS = $(SPRNG_HOME)/lib -lsprng -lstdc++ -lm
    SPRNG_INC = $(SPRNG_HOME)/include
  endif

endif

export FC HOST MACHINE LIBS MPIEXEC
