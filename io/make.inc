#Define message

MESSAGE = "Compiling IO modules"

#Define subdirectories for common code

#COMMONDIRS = ../interpolation/slatec/

#SUBDIRS = $(COMMONDIRS)/1d_spline \
#          $(COMMONDIRS)/2d_spline \
#          $(COMMONDIRS)/3d_spline \
#          $(COMMONDIRS)/spline_common \
#          $(COMMONDIRS)/xerror

#Define dependencies

graphics_io_mod.o: xdraw_io_mod.o