# Define message

MESSAGE = "Compiling f90 driver"

# Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/grid \
          $(COMMONDIR)/newtonGMRES

# Module order

MODS = setup_mod.F driver_mod.F

