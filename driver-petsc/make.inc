#Define message

MESSAGE = "Compiling PETSC driver"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/grid \
          $(COMMONDIR)/newtonGMRES

# Petsc include

include ${PETSC_DIR}/bmake/common/base
