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

#MODS = nksetup_mod.F ts_setup_mod.F ts_mod.F

#Define dependencies

ts_setup_mod.o: nksetup_mod.o
ts_mod.o: ts_setup_mod.o
