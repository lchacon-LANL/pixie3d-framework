#Define message

MESSAGE = "Compiling setup modules"

#Define subdirectories for common code

COMMONDIRS = ..

SUBDIRS = $(COMMONDIRS)/io \
          $(COMMONDIRS)/grid

#Define dependencies

iosetup_mod.o: setup_mod.o
