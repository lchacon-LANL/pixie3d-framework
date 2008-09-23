#Define message

MESSAGE = "Compiling driver setup modules"

#Define subdirectories for common code

COMMONDIRS = ..

SUBDIRS = $(COMMONDIRS)/io \
          $(COMMONDIRS)/grid\
          $(COMMONDIRS)/var_setup

#Define dependencies

