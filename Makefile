
# Read application specific variables

-include make.inc

# Read local configuration

ifndef ARCH
-include $(COMMONDIR)/make/make.mach.inc
endif

ifndef MODFLAG 
-include $(COMMONDIR)/make/make.comp.inc
endif

ifndef CONTRIBLIBS
-include $(COMMONDIR)/make/make.lib.inc
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

#FRMWRK_REL1=3
#FRMWRK_REL2=1
#override CPPFLAGS += $(PREPROC)FRMWRK_REL1=$(FRMWRK_REL1) $(PREPROC)FRMWRK_REL2=$(FRMWRK_REL2)

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
	-@rm -f *.o *.mod *.a

distclean: clean
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Main setup targets

setup: contrib_setup
	-@for subdir in `find . -name "make.inc" -exec dirname {} \;` ; do \
		rm $$subdir/makefile 2>/dev/null ; \
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
	$(MAKE) --no-print-directory -e -C contrib/lsode lib
	$(MAKE) --no-print-directory -e -C contrib/slatec lib
ifeq ($(ARPACK),t)
	$(MAKE) --no-print-directory -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack lib
ifdef BOPT
	$(MAKE) --no-print-directory -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack plib
endif
endif
	$(MAKE) --no-print-directory -e -C contrib/fpa/src lib
	$(MAKE) --no-print-directory -e -C contrib/sdc/src lib
	$(MAKE) --no-print-directory -e -C contrib/rng/src lib
	$(MAKE) --no-print-directory -e -C contrib/btridiag lib
ifeq ($(PIT),t)
	$(MAKE) --no-print-directory -e -C contrib/parareal all	
endif
ifeq ($(PTRID),t)
ifdef BOPT
	$(MAKE) --no-print-directory -e -C contrib/ptridiag lib
endif
endif

contrib_clean: ;
	@$(MAKE) -e -C contrib/lsode clean
	@$(MAKE) -e -C contrib/slatec distclean
ifeq ($(ARPACK),t)
	@$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack clean
endif
	@$(MAKE) -e -C contrib/fpa/src distclean
	@$(MAKE) -e -C contrib/sdc/src distclean
	@$(MAKE) -e -C contrib/parareal clean	
	@$(MAKE) -e -C contrib/rng/src clean
	@$(MAKE) -e -C contrib/btridiag clean
	@$(MAKE) -e -C contrib/ptridiag clean

contrib_setup: ;
	-@rm -r contrib
	-@tar xzf common_contrib.tgz

contrib_pack: ;
	-@rm -f common_contrib.tgz > /dev/null
	-@tar czf common_contrib.tgz contrib

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

%.o : %.C
	@echo 'Compiling' $@
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $<


