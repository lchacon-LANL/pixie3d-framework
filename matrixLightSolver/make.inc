#Define message

MESSAGE = "Compiling matrix-light solver modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math \
	  $(COMMONDIR)/grid
