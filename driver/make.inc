#Define message

MESSAGE = "Compiling f90 driver"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/grid \
          $(COMMONDIR)/newtonGMRES
