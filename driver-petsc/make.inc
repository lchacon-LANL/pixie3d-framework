#Define message

MESSAGE = "Compiling PETSC driver"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math $(COMMONDIR)/interpolation/1d_int $(COMMONDIR)/grid

# Petsc include
ifneq ($(HOSTNAME),dlg-hp.ornl.gov)
	include ${PETSC_DIR}/bmake/common/base
endif

# Module order

#MODS := setup_mod.F driver_mod.F
