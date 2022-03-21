FLAGS = -std=c18 -Wall -Wextra
ADDITIONAL_FLAGS =
HPC_FLAGS = -O3 -march=native

# OpenMP or MPI?
ifeq ($(code), mpi_A2A)
	SRC = main_mpi_A2A.c  kdtree_mpi_A2A.c kdtree_shared.c
	COMPILER = mpicc
	OUT = tree_mpi_A2A.x

else ifeq ($(code), openmp)
	SRC = main_openmp.c   kdtree_shared.c
	ADDITIONAL_FLAGS +=  -fopenmp
	COMPILER = gcc
	OUT = tree_openmp.x  

else ifeq ($(code), mpi)
	SRC = main_mpi.c  kdtree_mpi.c kdtree_shared.c
	COMPILER = mpicc
	OUT = tree_mpi.x
else
	SRC = main_serial.c kdtree_shared.c
	COMPILER = gcc
	OUT = tree_serial.x
endif


# Floating point or double precision?
ifeq ($(precision), double)
	ADDITIONAL_FLAGS += -D DOUBLE_PRECISION
endif

# loading data from file?
ifeq ($(data), from_file)
	ADDITIONAL_FLAGS += -D LOAD_DATA_FROM_FILE
endif

# version? Move or Store?
ifeq ($(version), store)
	ADDITIONAL_FLAGS += -D STORE
endif

default: $(SRC)
	$(COMPILER) $(HPC_FLAGS)  $(FLAGS) $(ADDITIONAL_FLAGS) -o $(OUT) $(SRC) -I.

debug_points: $(SRC)
	$(COMPILER) $(HPC_FLAGS)  $(FLAGS) $(ADDITIONAL_FLAGS) -D DEBUG_POINTS -o $(OUT) $(SRC) -I.

debug_tree: $(SRC)
	$(COMPILER) $(HPC_FLAGS)  $(FLAGS) $(ADDITIONAL_FLAGS) -D DEBUG_TREE -o $(OUT) $(SRC) -I.

print: $(SRC)
	$(COMPILER) $(HPC_FLAGS)  $(FLAGS) $(ADDITIONAL_FLAGS) -D PRINT_TREE  -o $(OUT) $(SRC) -I.

.PHONY: clean

clean:
	rm -f *.x

