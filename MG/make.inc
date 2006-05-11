#Define message

MESSAGE = "Compiling MG modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/matrixLightSolver
