#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/grid \
          $(COMMONDIR)/newtonGMRES

# Petsc include

include ${PETSC_DIR}/bmake/${PETSC_ARCH}/base
