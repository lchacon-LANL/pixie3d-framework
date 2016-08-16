# Define message

MESSAGE = "Compiling program driver"

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
          $(COMMONDIR)/eigen\
          $(COMMONDIR)/OMP

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

#Setup defs.

#LNK_FILES = ../driver-petsc/driver_mod.F 

#Dependencies
