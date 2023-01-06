#Library name

LIBNAME = libio.a

#Define message

MESSAGE = "Compiling IO modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS =

#Define dependencies

graphics_io_mod.o: io_mod.o
xdraw_io_mod.o: graphics_io_mod.o
