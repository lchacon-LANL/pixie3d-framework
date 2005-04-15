#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/matrixLightSolver-2d \
	  $(COMMONDIR)/newtonGMRES \
	  $(COMMONDIR)/fft \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/interpolation/slatec/1d_spline \
	  $(COMMONDIR)/interpolation/slatec/2d_spline \
	  $(COMMONDIR)/interpolation/slatec/xerror
