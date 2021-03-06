#Define message

MESSAGE = "Compiling plotting modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/matrixLightSolver\
          $(COMMONDIR)/MG\
          $(COMMONDIR)/OMP\
          $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

#Define dependencies

graphics_mod.o: xdraw_graphics_mod.o HDF5_graphics_mod.o ADIOS_graphics_mod.o ADIOS2_graphics_mod.o

3dplot.o: acc_test_mod.o

#Setup defs.

#LNK_FILES = ../driver-petsc/setup_mod.F
