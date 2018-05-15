#Library name

LIBNAME = libomp.a

#Define message

MESSAGE = "Compiling OMP modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/io \
	  $(COMMONDIR)/interpolation/1d_int

#Define dependencies

#xdraw_io_mod.o: graphics_io_mod.o
