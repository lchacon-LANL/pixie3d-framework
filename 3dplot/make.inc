#Define message

MESSAGE = "Compiling plotting modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/mk

#Define dependencies

graphics_mod.o: xdraw_graphics_mod.o HDF5_graphics_mod.o ADIOS_graphics_mod.o

3dplot.o: acc_test_mod.o

#Setup defs.

#LNK_FILES = ../driver-petsc/setup_mod.F
