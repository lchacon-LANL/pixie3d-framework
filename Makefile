# Read application specific variables

-include make.inc

ifndef MODFLAG 
-include $(COMMONDIR)/make/make.comp.inc
endif

# GENERAL PURPOSE MAKEFILE

SRC  ?= $(wildcard *.[F,c,f] *.f90 *.F90)

MODS ?= $(wildcard *_mod.[F,f] *_mod.f90 *_mod.F90)

OBJS   := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.C,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.f90,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F90,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(patsubst %.f,%.o,$(patsubst %.F,%.f,$(patsubst %.f90,%.f,$(patsubst %.F90,%.f,$(MODS)))))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F] $(dir)/*_mod.f90 $(dir)/*_mod.F90))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c,C] $(dir)/*.f90 $(dir)/*.F90)))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.C,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC))\
                            $(patsubst %.f90,%.o,$(COMMON_SRC))\
                            $(patsubst %.F90,%.o,$(COMMON_SRC)))

LIBS :=

prefix = .

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define main targets

PWD = `pwd`

.PHONY: setup setup_lnk target lib message clean distclean common contrib contrib_clean $(SUBDIRS)

target: common message $(OBJMOD) $(OBJS) $(COMMON_OBJS)

message: ;
	@echo ''
#	@echo $(MESSAGE)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

#Cleaning targets

clean: ;
	-rm -f *.o *.mod *.a

distclean: clean
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Main setup targets

setup: contrib_setup
	-@for subdir in `find . -name "make.inc" -exec dirname {} \;` ; do \
		-rm $$subdir/makefile 2>/dev/null ; \
		ln -s -f $(PWD)/Makefile $$subdir/makefile 2>/dev/null ; \
		$(MAKE) -C $$subdir setup_lnk; done

setup_lnk: ;
	-@for file in $(LNK_FILES) ; do \
		ln -s $$file 2>/dev/null ; done

#Library setup

lib: common $(OBJMOD) $(OBJS) $(COMMON_OBJS)
ifdef LIBNAME
	-ar rs $(LIBNAME) $(OBJMOD) $(OBJS) 
endif
ifdef LIBNAME_COM
	-ar rs $(LIBNAME_COM) $(OBJMOD) $(OBJS) $(COMMON_OBJS)
endif

#Contributed software setup

contrib: ;
	$(MAKE) -e -C contrib/lsode lib
	$(MAKE) -e -C contrib/slatec lib
ifdef ARPACK
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack lib
ifdef BOPT
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack plib
endif
endif
ifdef FPA
	$(MAKE) -e -C contrib/fpa/src lib
endif

contrib_clean: ;
	$(MAKE) -e -C contrib/lsode clean
	$(MAKE) -e -C contrib/slatec distclean
ifdef ARPACK
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack clean
endif
ifdef FPA
	$(MAKE) -e -C contrib/fpa/src distclean
endif

contrib_setup: ;
	-tar xzf common_contrib.tgz
	$(MAKE) -e -C contrib/slatec setup

#Define dependencies

$(OBJS) : $(MODS) $(COMMOM_MODS)

#Define patterns

%.o : %.f
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<

%.o : %.f90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<

%.o : %.c
	@echo 'Compiling' $@
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<
