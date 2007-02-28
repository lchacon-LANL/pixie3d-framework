# Define message

MESSAGE = "Compiling f90 driver"

# Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/grid \
          $(COMMONDIR)/newtonGMRES

# Module order

MODS = setup_mod.F driver_mod.F

#Setup defs.

LNK_FILES = ../driver-petsc/setup_mod.F\
            ../driver-petsc/driver_mod.F