#Define message

MESSAGE = "Compiling VARIABLE SETUP modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid

#Define dependencies

var_io_mod.o: var_setup_mod.o
var_bc_mod.o: var_setup_mod.o
