GSL_DIR = /software/gsl-2.2.1-el6-x86_64+intel-15.0
GSLPATH = /software/gsl-2.2.1-el6-x86_64+intel-15.0/lib
LIBLAPACKPATH = /project/gavoth/ynhan/local/lapack-3.5.0/
LAPACKEPATH = /project/gavoth/ynhan/local/lapack-3.5.0/lapacke/include
CBLASPATH = /project/gavoth/ynhan/local/CBLAS/include
LIBCBLAS = /project/gavoth/ynhan/local/CBLAS/lib/

WARN_FLAGS = -Wall

OPT = -O2 -std=c++11 $(WARN_FLAGS)
MKL_OPT = -O2 -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -fopenmp -std=c++11 $(WARN_FLAGS)

LIBS         = -lm -lgsl -mkl -lblas -llapacke -llapack -lcblas -lgfortran
LDFLAGS      = -L$(GSLPATH) -L$(LIBLAPACKPATH) -L$(LIBCBLAS)
CFLAGS	     = -I$(LAPACKEPATH) -I$(CBLASPATH)

CC           = icpc $(OPT)

test.x: FitGLE.o main.o Frame.o
	$(CC) $(LDFLAGS) -o FitGLE.x main.o FitGLE.o Frame.o $(LIBS)

main.o: main.cpp
	$(CC) -c main.cpp

FitGLE.o: FitGLE.cpp
	$(CC) $(CFLAGS) -c FitGLE.cpp

Frame.o: Frame.cpp
	$(CC) -c Frame.cpp


