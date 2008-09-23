#Define message

MESSAGE = "Compiling setup modules"

#Define subdirectories for common code

COMMONDIRS = ..

SUBDIRS = $(COMMONDIRS)/io \
          $(COMMONDIRS)/grid

#Define dependencies

var_io_mod.o: var_setup_mod.o
