#Define library name

LIBNAME = libslatec.a

#Define message

MESSAGE = "Compiling SLATEC interpolation libraries"

#Define subdirectories for common code

COMMONDIRS = .

SUBDIRS = $(COMMONDIRS)/1d_spline \
          $(COMMONDIRS)/2d_spline \
	  $(COMMONDIRS)/3d_spline \
	  $(COMMONDIRS)/spline_common \
	  $(COMMONDIRS)/xerror
