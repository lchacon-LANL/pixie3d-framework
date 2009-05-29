# Default

HOST?=$(HOSTNAME)

# Machine configuration

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
      HDF5 = t
    endif
    ifeq ($(FC),gfortran)
      HDF5_HOME =/usr/local/hdf5-1.6.7/gfortran/serial
      HDF5 = t
    endif
  endif
  HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf-3

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),quevedo.ornl.gov)
  VMEC = t
  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-2.3.3

#    PETSC_ARCH = linux64_absoft_
#    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/parallel

    FC=gfortran
    PETSC_ARCH = linux64_openmpi
    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
  else
    LIBS += -L/usr/lib64 -llapack -lblas
    HDF5 = f
    ifeq ($(FC),f90)
      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
      HDF5 = t
    endif
    ifeq ($(FC),gfortran)
      HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/serial
      HDF5 = t
    endif
  endif
  HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf-3

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),gongora.ornl.gov)
  FC=gfortran

  LIBS += -L/usr/lib -llapack -lblas

  PETSC_DIR =$(HOME)/lib/petsc-2.3.3
  PETSC_ARCH = linux_openmpi

  HDF5 = t
  HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
  HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm

  NETCDF=t
  NETCDF_LIBS = -lnetcdff -lnetcdf
  NETCDF_INC  = -I/usr/include/netcdf-3

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),kiwi)
  FC=g95

#  PETSC_DIR =$(HOME)/lib/petsc-2.3.3
#  PETSC_ARCH = linux_openmpi

  HDF5 = f
#   HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
#   HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm

#  NETCDF=t
#  NETCDF_LIBS = -lnetcdff -lnetcdf
#  NETCDF_INC  = -I/usr/include/netcdf-3

#  ARPACK_LIBS = -L../contrib/arpack/ -larpack_$(FC)
endif

ifeq ($(HOST),nip.lanl.gov)
  LIBS += -llapack -lblas -lg2c 
endif

ifeq ($(HOST),ra4)
  PETSC_DIR =/ricercatori/ft/petsc-2.2.0
  HDF5_HOME =
  CPPFLAGS += $(PREPROC)RFX
endif

ifeq ($(HOST),ra22)
  PETSC_ARCH = linux_openmpi
  PETSC_DIR =/ricercatori/ft/petsc-2.3.3
  HDF5_HOME =
  CPPFLAGS += $(PREPROC)RFX
  FC        = pgf95
  LIBS     += -llapack -lblas
endif

ifeq ($(HOST),cayenne1)
#  MPI_HOME  =/packages/mpi/mpich2-1.0.5p4-gcc-4.1-pgi-7.0-debug
  OPT       = g
  FC        = pgf95
  PETSC_C   = t
  BOPT      = O
#  FFLAGS    += -i8
endif

ifeq ($(HOST),gongora.lanl.gov)
  PETSC_DIR =/usr/local/petsc-2.2.0
  HDF5_HOME =/usr/local/hdf5/parallel/f90_
#  MPI_HOME  =/usr/local/mpich-1.2.5.2/f90_
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

ifeq ($(HOST),smokeylogin1)
  FC = pgf95
#  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS = $(HDF5_FLIB)
  HDF5 = t
  BOPT = t
endif