# Read application specific variables

-include make.local

#Default compiler

FC ?= mpif90
PREPROC ?= -D
MODFLAG ?= -I
ADDMODFLAG  ?= -I
VERBOSE     ?= -v
CPPFLAGS    ?= $(PREPROC)gfortran
FFLAGS      ?= -fno-automatic -O3

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

#Define targets

PWD = `pwd`

.PHONY: setup setup_lnk target lib message clean distclean common contrib contrib_clean $(SUBDIRS)

lib: common $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	-ar rs $(LIBNAME) $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	-mv $(LIBNAME) $(LIB_PATH)
	-cp *.mod *.MOD $(LIB_PATH) 2>/dev/null || true

target: common message $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	-cp *.mod *.MOD $(LIB_PATH) 2>/dev/null || true

message: ;
	@echo ''
#	@echo $(MESSAGE)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target LIB_PATH=../$(LIB_PATH)

clean: ;
	-@rm -f *.o *.mod *.MOD *.a lib/*.a lib/*.mod lib/*.MOD || true

distclean: clean
	-@rm lib/*.a lib/*.mod lib/*.MOD 2>/dev/null || true
	-@for subdir in $(SUBDIRS) ; do \
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

%.o : %.f90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<

%.o : %.c
	@echo 'Compiling' $@
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<
