#Define library name

LIBNAME = libslatec.a

#Define message

MESSAGE = "Compiling SLATEC interpolation libraries"

#Define subdirectories for common code

COMMONDIR = ../..

SUBDIRS = ./1d_spline \
          ./2d_spline \
	  ./3d_spline \
	  ./spline_common \
	  ./xerror
