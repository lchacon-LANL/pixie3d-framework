#Define message

MESSAGE = "Compiling setup modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/io \
          $(COMMONDIR)/grid

#Define dependencies

var_io_mod.o: var_setup_mod.o
