# Default

HOST?=$(HOSTNAME)

# Machine configuration

ifeq ($(HOST),quevedo.ornl.gov)
  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-3.0.0-p7

    FC=gfortran
    ifeq ($(OPT),O)
      PETSC_ARCH = linux64_openmpi
    else
      PETSC_ARCH = linux64_openmpi_g
    endif

    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
    HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
    HDF5_MOD  = $(HDF5_HOME)/lib
    HDF5_INC  =-I$(HDF5_HOME)/include
  else
    LIBS = -L/usr/lib64 -llapack -lblas
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

#      HDF5_HOME =/usr/lib64
#      HDF5_LIBS =-L$(HDF5_HOME) -lhdf5_fortran -lhdf5 -lz -lm
#      HDF5_MOD  =/usr/lib64/gfortran/modules
#      HDF5_INC  =-I/usr/include
      HDF5 = t
    endif
  endif

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf

  VMEC = t
#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),gongora.ornl.gov)
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
  VMEC = t
  ifdef BOPT
    PETSC_DIR =/usr/local/petsc-2.3.3

#    PETSC_ARCH = linux64_absoft_
#    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/parallel

    FC=gfortran
    PETSC_ARCH = linux64_openmpi
    HDF5_HOME =/usr/local/hdf5-1.6.7/gfortran/parallel

#    MPI_HOME = /usr
  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f
    ifeq ($(FC),f90)
      HDF5_HOME =/usr/local/hdf5-1.6.7/absoft_/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include
      HDF5 = t
    endif
    ifeq ($(FC),gfortran)
      HDF5_HOME =/usr/local/hdf5-1.6.7/gfortran/serial
      HDF5_LIBS =-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
      HDF5_MOD  = $(HDF5_HOME)/lib
      HDF5_INC  =-I$(HDF5_HOME)/include
      HDF5 = t
    endif
  endif

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf-3

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),ra22)
  PETSC_ARCH = linux_openmpi
  PETSC_DIR =/ricercatori/ft/petsc-2.3.3
  HDF5_HOME =
  CPPFLAGS += $(PREPROC)RFX
  FC        = pgf95
  LIBS     += -llapack -lblas
endif

ifeq ($(HOST),bassi)
  OPTIMIZATION = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
  DEBUG        = -g
  PROFILE      = -P
  STATIC       = -s
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
#  FFLAGS       = -qmaxmem=-1 -bmaxstack:0x22000000 -qsave -blpdata -qsmallstack -qfixed=72
  FFLAGS       = -qfixed=72

  PREPROC      = -WF,-D
  CPPFLAGS    += $(PREPROC)xlf

  BOPT = t
endif

ifeq ($(HOST),franklin)
  FC = pgf95
  HDF5_HOME = $(HDF5_PAR_DIR)
#  HDF5_LIBS:= $(HDF5)
  override HDF5 = t
  BOPT = t
endif

ifeq ($(findstring smoky,$(HOST)),smoky)
  FC = pgf95
#  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5 = t
  BOPT = t
endif

ifeq ($(findstring jaguar,$(HOST)),jaguar)
  ifeq ($(PE_ENV),GNU) 
    FC = gfortran
  endif
  ifeq ($(PE_ENV),PGI)
    FC = pgf95
  endif

  HDF5 = t
#  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS = $(HDF5_FLIB)

  BOPT = t
  PETSC_DIR ?= /opt/petsc/2.3.3/real/${PE_ENV}/${XTPE_COMPILE_TARGET}
  PETSC_ARCH ?= cray-xt
endif
