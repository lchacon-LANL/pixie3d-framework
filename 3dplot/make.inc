#Define message

MESSAGE = "Compiling plotting modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math  $(COMMONDIR)/interpolation/1d_int $(COMMONDIR)/grid 

#Define dependencies

xdraw_graphics_mod.o: setup_mod.o

HDF5_graphics_mod.o: setup_mod.o plot_mod.o

graphics_mod.o: xdraw_graphics_mod.o HDF5_graphics_mod.o

