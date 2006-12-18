#Define message

MESSAGE = "Compiling MG modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/interpolation/1d_int \
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/matrixLightSolver
