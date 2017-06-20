#Define message

MESSAGE = "Compiling driver setup modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io \
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid\
          $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/newtonGMRES \
          $(COMMONDIR)/eigen \
          $(COMMONDIR)/OMP

# Module order

#Define dependencies

ts_setup_mod.o: ts_nk_setup_mod.o
ts_base_mod.o: ts_setup_mod.o
ts_io_mod.o: ts_base_mod.o
ts_if_mod.o: ts_io_mod.o
proj_int_mod.o: ts_if_mod.o
ts_mod.o: proj_int_mod.o
