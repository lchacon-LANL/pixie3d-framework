#Define message

MESSAGE = "Compiling orbit modules"

#Library name

LIBNAME = liborbit.a

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/grid\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/interpolation/slatec

#Define dependencies

orbit_mod.o: spline_mod.o
