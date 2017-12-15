#Define message

MESSAGE = "Compiling orbit modules"

#Library name

LIBNAME = liborbit.a

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
	  $(COMMONDIR)/OMP\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup

#Define dependencies

orbit_mod.o: lyapn_mod.o
par_int_mod.o: orbit_mod.o poincare_io_mod.o
