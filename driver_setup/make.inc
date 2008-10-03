#Define message

MESSAGE = "Compiling driver setup modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/io \
          $(COMMONDIR)/grid\
          $(COMMONDIR)/var_setup

#Define dependencies

