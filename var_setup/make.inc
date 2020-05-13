#Define message

MESSAGE = "Compiling VARIABLE SETUP modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid

#Define dependencies

var_bc_mod.o: var_setup_mod.o
var_mod.o: var_bc_mod.o
ADIOS_io_mod.o: var_mod.o
ADIOS2_io_mod.o: var_mod.o
var_io_mod.o: ADIOS_io_mod.o ADIOS2_io_mod.o
