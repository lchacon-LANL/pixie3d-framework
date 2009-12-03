#Define message

MESSAGE = "Compiling f90 NK modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math $(COMMONDIR)/grid 

#Dependencies

nk_mod.o: fixed_point_accelerator.o
fixed_point_accelerator.o: fpa_kinds.o
fpa_kinds.o: f90_assert.o

#Specific code

CPPFLAGS += -DNDEBUG
