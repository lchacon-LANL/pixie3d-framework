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

#Module search path

MODDIRS = $(GMODPATH) $(MODPATH). $(patsubst $(COMMONDIR)%,$(ADDMODPATH)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

.PHONY: message clean distclean common $(SUBDIRS)

target: common message $(OBJMOD) $(OBJS) $(COMMON_OBJS)

message: ;
	@echo ''
#	@echo $(MESSAGE)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

clean:
	-rm -f *.o *.mod

distclean: clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS) : $(MODS) $(COMMOM_MODS)

#Define patterns

%.o : %.f
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<

#%.o : %.c
#	@echo 'Compiling' $@
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<
