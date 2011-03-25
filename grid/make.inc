#Define message

MESSAGE = "Compiling grid modules"

#Library name

LIBNAME = libgrid.a

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math $(COMMONDIR)/io $(COMMONDIR)/interpolation/1d_int

#Define dependencies

grid_pack_mod.o: grid_petsc_mod.o 

grid_def_mod.o: grid_pack_mod.o

grid_operations_mod.o: grid_def_mod.o

grid_anal_map_mod.o: grid_def_mod.o grid_operations_mod.o

grid_mpi_mod.o: grid_def_mod.o

grid_mg_mod.o: grid_mpi_mod.o grid_operations_mod.o

grid_bc_mod.o: grid_mg_mod.o grid_operations_mod.o

grid_metric_mod.o: grid_bc_mod.o grid_anal_map_mod.o grid_debug_mod.o

grid_create_mod.o: grid_metric_mod.o

grid_diffops_mod.o: grid_create_mod.o

grid_mod.o: grid_create_mod.o grid_operations_mod.o grid_diffops_mod.o
