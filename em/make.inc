#Define message

MESSAGE = "Compiling EM modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/OMP\
	  $(COMMONDIR)/pic/2d
	  
#Define dependencies

regularization_aa_mod.o: aamem_utilities_mod.o
gm_utils_mod.o: regularization_aa_mod.o
gs_kmeans_mod.o: gm_utils_mod.o
kmeans_mod.o: gs_kmeans_mod.o
gm_init_mod.o: kmeans_mod.o
EM_utils_mod.o: gm_init_mod.o
EM_mod.o: EM_utils_mod.o
AMEM_mod.o: EM_utils_mod.o

