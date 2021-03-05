#Define message

MESSAGE = "Compiling 2d PIC modules"

#Define subdirectories for common code

COMMONDIR = ../..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io \
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid\
          $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
          $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
          $(COMMONDIR)/OMP\
          $(COMMONDIR)/orbit

#Module order

#Whether canonical momenta diag are compiled
PIC_CM_DIAG?=t

ifeq ($(PIC_CM_DIAG),t) 
  override CPPFLAGS += $(PREPROC)PIC_CM_DIAG
endif

#Obsolete
PIC_CUDA=f
USE_TEX=f
GOTOBLAS = f
#intel math lib
LIBIMF = t

#for c++
ifeq ($(CC),g++)
    Cpp_OPT=-lstdc++
endif

#for pic cuda 
ifeq ($(PIC_CUDA),t) 
    override CPPFLAGS += $(PREPROC)_CUDA
    ifeq ($(USE_TEX),t) 
       override CPPFLAGS += $(PREPROC)_USE_TEX
    endif    
    SRC_cu = $(wildcard *.cu)
    OBJ_cu = $(SRC_cu:.cu=.o)
endif

#Define dependencies

Fcc_mod.o: pic_defs.h
pic_base_mod.o: Fcc_mod.o
pic_mod.o: pic_base_mod.o

ifeq ($(PIC_CUDA),t) 
pic.o : pic.cu pic_kernel.cuh
	@echo 'Compiling' $@
	$(NVCC) -c ${gpu_arch} $(CUDA_INCLUDE) $(NVCCFLAGS) $(CPPFLAGS) $<
Fcc_mod.o   : pic.o
endif
