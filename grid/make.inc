#Define message

MESSAGE = "Compiling grid modules"

#Define subdirectories for common code

COMMONDIR = ..

SUBDIRS = $(COMMONDIR)/math $(COMMONDIR)/interpolation/1d_int

#Define dependencies

grid_debug_mod.o: xdraw_io_mod.o

grid_operations_mod.o: grid_def_mod.o

grid_anal_map_mod.o: grid_def_mod.o grid_operations_mod.o

grid_mpi_mod.o: grid_def_mod.o

grid_mg_mod.o: grid_mpi_mod.o grid_operations_mod.o

grid_bc_mod.o: grid_mg_mod.o grid_operations_mod.o

grid_metric_mod.o: grid_bc_mod.o grid_anal_map_mod.o grid_debug_mod.o

grid_create_mod.o: grid_metric_mod.o

grid_mod.o: grid_create_mod.o grid_operations_mod.o