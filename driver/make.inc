# Define message

MESSAGE = "Compiling f90 driver"

# Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/io \
	  $(COMMONDIR)/grid \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/var_setup \
	  $(COMMONDIR)/driver_setup \
          $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/matrixLightSolver\
          $(COMMONDIR)/MG\
          $(COMMONDIR)/orbit\
          $(COMMONDIR)/mk\
          $(COMMONDIR)/eigen

# Module order

MODS = newtongm_mod.F

#Setup defs.

#LNK_FILES = ../driver-petsc/driver_mod.F 

#Dependencies
