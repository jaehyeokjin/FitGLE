LIBLAPACKPATH = /project/gavoth/ynhan/local/lapack-3.5.0/
LAPACKEPATH = /project/gavoth/ynhan/local/lapack-3.5.0/lapacke/include
OPT = -O3 -Wall

CFLAGS = $(OPT) -I$(LAPACKEPATH)
LDFLAGS = $(OPT) -L$(LIBLAPACKPATH)
LIBS = -llapacke -llapack -lblas -lgfortran
CC = gcc

test.x : dgel.o lapacke_example_aux.o
	$(CC) $(LDFLAGS) -o $@ dgel.o lapacke_example_aux.o  $(LIBS)

dgel.o : dgel.c
	$(CC) $(CFLAGS) -c dgel.c

lapacke_example_aux.o : lapacke_example_aux.c
	$(CC) $(CFLAGS) -c lapacke_example_aux.c



