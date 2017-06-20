#Define message

MESSAGE = "Compiling matrix-light solver modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/io \
	  $(COMMONDIR)/math \
	  $(COMMONDIR)/grid \
	  $(COMMONDIR)/interpolation/1d_int
