GSL_DIR = /software/gsl-2.2.1-el6-x86_64+intel-15.0
GSLPATH = $(GSL_DIR)/lib

WARN_FLAGS = -Wall

OPT = -O2 -std=c++11 $(WARN_FLAGS)
MKL_OPT = -O2 -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -fopenmp -std=c++11 $(WARN_FLAGS)

LIBS         = -lm -L$(GSLPATH) -lgsl -mkl
LDFLAGS      = $(OPT) 
CFLAGS	     = $(OPT)

CC           = icc

test.x: gsl_test.o
	$(CC) $(LDFLAGS) -o test.x gsl_test.o $(LIBS)

gsl_test.o: gsl_test.cpp
	$(CC) $(CFLAGS) -c gsl_test.cpp


