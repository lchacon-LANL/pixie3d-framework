# Read application specific variables

include make.inc

# GENERAL PURPOSE MAKEFILE

SRC  := $(wildcard *.[F,f,c])

MODS := $(wildcard *_mod.[F,f])

OBJS  := $(filter %.o,$(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(filter %.o, $(patsubst %.f,%.o,$(MODS))\
                        $(patsubst %.F,%.o,$(MODS)))


COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F]))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c])))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC)))

LIBS :=

prefix = .

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  MODPATH = -p
  MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH)$(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  MODPATH = -M
  MODDIRS = $(MODPATH) . $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  MODPATH =
  MODDIRS = $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

#Define targets

.PHONY: clean recclean common $(SUBDIRS)

target: common $(OBJMOD) $(OBJS) $(COMMON_OBJS)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

clean:
	-rm -f *.o *.mod

recclean: clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS) : $(MODS) $(COMMOM_MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	$(FC) -c $(MODDIRS) $(FFLAGS) $<
