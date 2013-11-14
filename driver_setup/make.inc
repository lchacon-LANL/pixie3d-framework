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
          $(COMMONDIR)/eigen

# Module order

#Define dependencies

ts_setup_mod.o: nksetup_mod.o
ts_base_mod.o: ts_setup_mod.o
proj_int_mod.o: ts_base_mod.o
ts_mod.o: proj_int_mod.o
