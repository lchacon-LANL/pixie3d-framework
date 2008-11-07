#Define message

MESSAGE = "Compiling f90 eigen modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS =$(COMMONDIR)/math $(COMMONDIR)/grid $(COMMONDIR)/newtonGMRES
