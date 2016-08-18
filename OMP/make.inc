#Library name

LIBNAME = libomp.a

#Define message

MESSAGE = "Compiling OMP modules"

#Define subdirectories for common code

COMMONDIR = $(PWD)/common

SUBDIRS = $(COMMONDIR)/grid

#Define dependencies

#xdraw_io_mod.o: graphics_io_mod.o
