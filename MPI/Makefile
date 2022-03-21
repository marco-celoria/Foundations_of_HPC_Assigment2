CC     = mpicc
CFLAGS = -Wall -O3 -march=native -std=c11
LDFLAGS = 

SRC = KDtree_MPI_last2.c
EXE = $(SRC:.c=)
OUT = $(SRC:.c=.out)

.PHONY:	default clean

default: $(EXE)

%:	%.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm $(EXE)
