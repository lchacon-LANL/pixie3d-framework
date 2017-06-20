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
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/OMP

#Define dependencies

spline_mod.o: lyapn_mod.o
orbit_mod.o: spline_mod.o
par_int_mod.o: orbit_mod.o poincare_io_mod.o
