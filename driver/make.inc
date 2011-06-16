# Define message

MESSAGE = "Compiling f90 driver"

# Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/io \
	  $(COMMONDIR)/grid \
	  $(COMMONDIR)/var_setup \
	  $(COMMONDIR)/driver_setup \
          $(COMMONDIR)/newtonGMRES\
          $(COMMONDIR)/MG\
          $(COMMONDIR)/mk\
          $(COMMONDIR)/eigen

# Module order

MODS = newtongm_mod.F

#Setup defs.

#LNK_FILES = ../driver-petsc/driver_mod.F 

#Dependencies
