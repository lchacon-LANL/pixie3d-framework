#Define message

MESSAGE = "Compiling plotting modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup

#Define dependencies

#xdraw_graphics_mod.o: setup_mod.o

#HDF5_graphics_mod.o: setup_mod.o plot_mod.o

graphics_mod.o: xdraw_graphics_mod.o HDF5_graphics_mod.o

#Setup defs.

#LNK_FILES = ../driver-petsc/setup_mod.F
