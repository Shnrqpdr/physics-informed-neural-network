CC = gcc
CFLAGS = -lm
OBJ = schrodinger-solver.o

ODIR = solver

prg: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o