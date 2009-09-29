# Read application specific variables

-include make.inc

ifndef MODFLAG 
-include $(COMMONDIR)/make/make.comp.inc
endif

# GENERAL PURPOSE MAKEFILE

SRC  ?= $(wildcard *.[F,f,c])

MODS ?= $(wildcard *_mod.[F,f])

OBJS  := $(filter %.o,$(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.C,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD :=  $(patsubst %.f,%.o,$(patsubst %.F,%.f,$(MODS)))


COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F]))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c])))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.C,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC)))

LIBS :=

prefix = .

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

PWD = `pwd`

.PHONY: setup setup_lnk target lib message clean distclean common $(SUBDIRS)

target: common message $(OBJMOD) $(OBJS) $(COMMON_OBJS)

message: ;
	@echo ''
#	@echo $(MESSAGE)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

clean: ;
	-rm -f *.o *.mod *.a

distclean: clean
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

setup: ;
	-@for subdir in `find . -name "make.inc" -exec dirname {} \;` ; do \
		-rm $$subdir/makefile 2>/dev/null ; \
		ln -s -f $(PWD)/Makefile $$subdir/makefile 2>/dev/null ; \
		$(MAKE) -C $$subdir setup_lnk; done

setup_lnk: ;
	-@for file in $(LNK_FILES) ; do \
		ln -s $$file 2>/dev/null ; done

lib: common $(OBJMOD) $(OBJS) $(COMMON_OBJS)
ifdef LIBNAME
	-ar rs $(LIBNAME) $(OBJMOD) $(OBJS) 
endif
ifdef LIBNAME_COM
	-ar rs $(LIBNAME_COM) $(OBJMOD) $(OBJS) $(COMMON_OBJS)
endif

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
