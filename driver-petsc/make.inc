#Define message

MESSAGE = "Compiling PETSC driver"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math $(COMMONDIR)/interpolation/1d_int $(COMMONDIR)/grid

# Petsc include

include ${PETSC_DIR}/bmake/common/base

# Module order

MODS := setup_mod.F driver_mod.F
