# Define message

MESSAGE = "Compiling SAMRAI driver"

# Define subdirectories for common code

COMMONDIR = ..

#SUBDIRS = $(COMMONDIR)/math \
#	  $(COMMONDIR)/interpolation/1d_int \
#	  $(COMMONDIR)/io \
#	  $(COMMONDIR)/grid \
#	  $(COMMONDIR)/var_setup \
#	  $(COMMONDIR)/driver_setup \
#          $(COMMONDIR)/newtonGMRES\
#         $(COMMONDIR)/eigen

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/io \
	  $(COMMONDIR)/grid \
	  $(COMMONDIR)/var_setup \
	  $(COMMONDIR)/driver_setup\
      $(COMMONDIR)/newtonGMRES \

# Module order

MODS = newtongm_mod.F
#MODS =  setup_mod.F driver_mod.F

#Setup defs.

#LNK_FILES = ../driver-petsc/driver_mod.F 

#Dependencies
