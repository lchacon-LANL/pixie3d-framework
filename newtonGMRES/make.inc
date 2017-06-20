#Define message

MESSAGE = "Compiling f90 NK modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/io $(COMMONDIR)/math $(COMMONDIR)/grid $(COMMONDIR)/interpolation/1d_int

#Dependencies

#Specific code
