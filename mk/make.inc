#Define message

MESSAGE = "Compiling Monge-Kantorovich modules"

#Library name

LIBNAME = libmk.a

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/grid\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/var_setup\
          $(COMMONDIR)/driver_setup\
          $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/matrixLightSolver\
          $(COMMONDIR)/MG\
          $(COMMONDIR)/orbit

#Define dependencies

mk_pc_mod.o: mk_mod.o
