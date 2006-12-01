#Define message

MESSAGE = "Compiling grid modules"

#Define subdirectories for common code

SUBDIRS = ../math

#Define dependencies

grid_def_mod.o: grid_structures_mod.o

bc_def_mod.o: grid_structures_mod.o

local_mpi_mod.o: bc_def_mod.o grid_def_mod.o

grid_metric_mod.o: grid_def_mod.o local_mpi_mod.o bc_def_mod.o

grid_operations_mod.o: grid_metric_mod.o

grid_mod.o: local_mpi_mod.o grid_operations_mod.o